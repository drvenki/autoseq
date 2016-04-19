import os


def normpath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def stripsuffix(thestring, suffix):
    """
    Strip the given suffix from the given string
    :param thestring: string from which to strip
    :param suffix: suffix to strip
    :return: string without the given suffix. If the given string
    doesn't end in suffix, the original string is returned
    """
    if thestring.endswith(suffix):
        return thestring[:-len(suffix)]
    return thestring
