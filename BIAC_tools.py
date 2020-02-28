#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Small tools for multiprocessing
"""

import multiprocessing, io, datetime, os, smtplib
from email.mime.text import MIMEText

from types import ModuleType, FunctionType
from gc import get_referents
import sys

mylogin = "jas297"

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType

class MyPool(multiprocessing.pool.Pool):
    #This whole class allows us to run sub processes within subprocesses (subprocesses cannot exceed max count though)
    Process = NoDaemonProcess


def isempty(object):
    if object is None:
        return True
    elif len(object) == 0:
        return True
    else:
        return False

def getsize(obj):
    #sum size of object & all members within that object.
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size

def send_mail(msg_txt,subject="Cluster message"):
    #Send mail with specified txt (and subject) to default address specified in global variable
    msg_content = io.StringIO()
    msg_content.write("Datetime : %s \n\n" % datetime.datetime.now())
    msg_content.write("JobID : %d \n\n" % os.getpid() )
    msg_content.write("Message : %s \n\n" % msg_txt)
    msg = MIMEText(msg_content.getvalue())
    msg_content.close()
    to_addr = "%s@duke.edu" % mylogin
    from_addr = "%s@duke.edu" % mylogin
    msg['Subject'] = subject
    msg['from'] = from_addr
    msg['to'] = to_addr
    s = smtplib.SMTP('smtpgw.duhs.duke.edu')
    s.sendmail(from_addr, [to_addr], msg.as_string())
    s.quit()
