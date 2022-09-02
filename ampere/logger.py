""" 
A module for controlling logging in ampere

This module defines all ampere's functionality regarding logging. It is 
intended to be imported in other parts of ampere and not accessed directly by
users.
"""



import sys
import logging
from datetime import datetime

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
    logger = logging.getLogger("ampere logger")
    logger.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))



class Logger(object):
    """
    A basic class to control logging. 

    Can be used as a mixin or attribute of another class, or called standalone.
    Should define methods to set up, modify, or terminal logging behaviour.
    """
    def setup_logging(self, verbose=False, to_file=True, to_terminal = True, logfile=None, logPath = None):
        if logfile is None: #beware! because of the below check for existing FileHandlers this will only get used the first time the logger is setup.
            logfile = "ampere_"+str(datetime.now()).replace(' ','_').replace(":","-")[:-7] #we santise the time stamp and remove the fractions of a second, since we are unlikely to run ampere in the same place more than once in under a second given how long it takes.
        if logPath is None:
            logPath = '.'
        level = logging.INFO
        logger = logging.getLogger("ampere logger")
        logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-8.8s]  %(message)s")
        self.logger = logger
        self.logger.propagate=False
        if verbose:
            level = logging.DEBUG

        self.logger.setLevel(level)
        if not len(logger.handlers): #check if any handlers have been defined before so we don't end up duplicating outputs
            if to_file and not any([type(h) is logging.FileHandler for h in logger.handlers]):
                fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, logfile))
                fileHandler.setFormatter(logFormatter)
                self.logger.addHandler(fileHandler)
            if to_terminal and not any([type(h) is logging.StreamHandler for h in logger.handlers]): #slightly problematic here if someone wants to define output to a different stream than stdout, but right now they would have to do that manually anyway
                handler = logging.StreamHandler(stream=sys.stdout)
                handler.setFormatter(logFormatter)
                self.logger.addHandler(handler)
        
        sys.excepthook = handle_exception












if __name__ == "__main__":
    Logger().setup_logging()
    raise RuntimeError("Test unhandled")
    print("This code didn't exit")
