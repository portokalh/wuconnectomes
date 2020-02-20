#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Small tools for multiprocessing
"""

import multiprocessing, io, datetime, os, smtplib
from email.mime.text import MIMEText
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


class MyPool(multiprocessing.pool.Pool):
    #This whole class allows us to run sub processes within subprocesses (subprocesses cannot exceed max count though)
    Process = NoDaemonProcess



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
