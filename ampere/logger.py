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
    Log uncaught exceptions and make sure they are correctly logged before Python exits.

    This function is called when an uncaught exception occurs and handles the exception by logging it with the logger
    "ampere logger". If the exception is a KeyboardInterrupt, the original exception hook is used to handle it.

    This function was taken from the following stack overflow answer: https://stackoverflow.com/a/16993115/16164384

    Parameters
    ---------
    exc_type : exception type
    Type of the exception.
    exc_value : exception instance
    Instance of the exception.
    exc_traceback : traceback
    Traceback of the exception.

    Generated with Chat-GPT
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
        """
        Set up logging for the application.

        This function configures the logger "ampere logger" to output log messages at the specified level. The level can be
        either logging.DEBUG or logging.INFO, depending on the value of the verbose parameter. By default, the level is
        logging.INFO.

        The logger can output log messages to a file, to a terminal, or both. The to_file and to_terminal parameters control
        where the log messages are output. If to_file is True, the log messages are output to a file. If to_terminal is
        True, the log messages are output to the terminal. If both are True, the log messages are output to both the file
        and the terminal.

        The file and terminal outputs use different formatters. The file formatter is "%(asctime)s [%(threadName)-12.12s]
        [%(levelname)-8.8s] %(message)s", while the terminal formatter is "%(levelname)-8.8s %(message)s".

        If the logfile parameter is not specified, the logger creates a log file with a name in the format "ampere_timestamp",
        where "timestamp" is the current date and time with the space replaced by an underscore and the colons replaced by
        dashes. If the logPath parameter is not specified, the logger uses the current directory as the path for the log file.

        This function also sets the exception hook to handle_exception to make sure that uncaught exceptions are correctly
        logged before Python exits.

        Parameters
        ----------
        verbose : bool, optional
        If True, the logger outputs log messages at the logging.DEBUG level. If False, the logger outputs log
        messages at the logging.INFO level. Default is False.
        to_file : bool, optional
        If True, the logger outputs log messages to a file. Default is True.
        to_terminal : bool, optional
        If True, the logger outputs log messages to the terminal. Default is True.
        logfile : str, optional
        The name of the log file. If not specified, the logger creates a log file with a name in the format "ampere_timestamp".
        logPath : str, optional
        The path to the log file. If not specified, the logger uses the current directory as the path.

        Generated with Chat-GPT
        """
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
