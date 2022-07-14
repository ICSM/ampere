""" 
A module for controlling logging in ampere

This module defines all ampere's functionality regarding logging. It is 
intended to be imported in other parts of ampere and not accessed directly by
users.
"""



import sys
import logging

logger = logging.getLogger("ampere logger")
logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-8.8s]  %(message)s")

def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Make sure that uncaught exceptions are correctly logged before python exits.

    This function was taken from the following stack overflow answer: https://stackoverflow.com/a/16993115/16164384
    """
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))



class Logger(object):
    """
    A basic class to control logging. 

    Can be used as a mixin or attribute of another class, or called standalone.
    Should define methods to set up, modify, or terminal logging behaviour.
    """
    def setup_logging(self, verbose=False, to_file=True, to_terminal = True, logfile="ampere", logPath = '.'):
        level = logging.INFO
        if verbose:
            level = logging.DEBUG

        logger.setLevel(level)
        if to_file:
            fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, logfile))
            fileHandler.setFormatter(logFormatter)
            logger.addHandler(fileHandler)
        if to_terminal:
            handler = logging.StreamHandler(stream=sys.stdout)
            handler.setFormatter(logFormatter)
            logger.addHandler(handler)
        
        sys.excepthook = handle_exception












if __name__ == "__main__":
    Logger().setup_logging()
    raise RuntimeError("Test unhandled")
    print("This code didn't exit")
