import json
import os
import datetime
import socket
import pymongo
import logging

import signal

from autoseq.util.path import normpath
from autoseq.util.status import Status


def remotelog(func):
    def run_logged_func(*args, **kwargs):
        return_value = None
        stats = None
        if 'sampledata' in kwargs:
            # add sample to logging server
            sampledata = kwargs['sampledata']
            mongolog = AutoseqLogging("~/.autoseqlogging.json", sampledata['SUBSTUDY'])
            sampledata[AutoseqLogging.STATUS_KEY] = Status.RUNNING
            sampledata[AutoseqLogging.START_TIME_KEY] = str(datetime.datetime.now())
            sampledata[AutoseqLogging.HOSTNAME_KEY] = socket.gethostname()
            sampledata[AutoseqLogging.LOGFILE_KEY] = None
            mongolog.add_run(sampledata)

            def capture_sigint(sig, frame):
                """
                Capture ctrl-c (or SIGINT sent in other ways).
                1. update remote log
                :param sig:
                :param frame:
                :return:
                """

            # signal.signal(signal.SIGINT, capture_sigint)
            try:
                return_value, stats = func(*args, **kwargs)
            except OSError:
                sampledata[AutoseqLogging.STATUS_KEY] = Status.CANCELLED
                mongolog.update_run(sampledata)
                logging.info("Autoseq was cancelled.")
                raise OSError

            if sampledata[AutoseqLogging.STATUS_KEY] == Status.RUNNING:
                sampledata[AutoseqLogging.STATUS_KEY] = Status.COMPLETED

            sampledata[AutoseqLogging.END_TIME_KEY] = str(datetime.datetime.now())
            mongolog.update_run(sampledata)
        else:
            raise ValueError("sampledata not supplied to logged function.")

        return return_value, stats

    return run_logged_func


class AutoseqLogging(object):
    START_TIME_KEY = "START_TIME"
    END_TIME_KEY = "END_TIME"
    SUBSTUDY_KEY = "SUBSTUDY"
    REPORTID_KEY = "REPORTID"
    STATUS_KEY = "STATUS"
    HOSTNAME_KEY = "HOSTNAME"
    LOGFILE_KEY = "LOGFILE"

    def __init__(self, credentials_file, collection_name):
        self.collection = None
        self.objid = None
        self.credentials_file = normpath(credentials_file)
        self.connect(collection_name)

    def connect(self, collection_name):
        """
        open a connection to our logging server and set jobs collection
        :param collection_name:
        :param credentials_file:
        """
        if os.path.exists(self.credentials_file):
            logging.info("Credentials file found, will log status to remote server")
            cred = json.load(open(self.credentials_file))
            uri = cred['uri'].format(user=cred['user'], pwd=cred['pwd'])
            client = pymongo.MongoClient(uri)
            db = client.get_default_database()
            self.collection = db[collection_name]
        else:
            logging.info("Credentials file {} not found. Can't log remotely.".format(self.credentials_file))

    def add_run(self, sampledata):
        """
        Add a run to the remote mongodb
        :param sampledata: dict to insert
        :rtype: None or ObjectID of inserted object (also set in self.objid)
        """
        if self.collection:
            logging.debug("Will insert {}".format(sampledata))
            try:
                self.objid = self.collection.insert(sampledata)
                logging.info("Inserted data with _id {}".format(self.objid))
            except pymongo.errors.ServerSelectionTimeoutError:
                logging.info("Count'n connect to mongodb server. Not logging.")
                self.objid = None
            return self.objid
        else:
            self.objid = None
            return None

    def update_run(self, sampledata):
        """
        :param sampledata: dict to update
        """
        if self.collection and self.objid:
            query = {'_id': self.objid}
            found_jobs = self.collection.find(query)
            job_count = found_jobs.count()
            if job_count == 0:  # objid not found
                logging.info("Object with _id {} not found. Could not update it. ".format(self.objid))
            elif job_count == 1:  # Update
                logging.debug("Will update _id {} with data {}".format(self.objid, sampledata))
                self.collection.update(query, {'$set': sampledata})
                logging.info("Updated {}".format(self.objid))

            if job_count > 1:  # Err out, error has occured when submitting
                raise ValueError("Combination of REPORTID and SUBSTUDY must be unique")
        else:
            pass


