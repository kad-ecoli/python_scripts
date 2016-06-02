#!/usr/bin/env python
# 2015-12-03 Chengxin Zhang
# read_mail.py 2
#   read the two most recent system mail
import os,sys
mail=open("/var/spool/mail/"+os.getenv("USER"),'rU').read()
if len(sys.argv)==1:
    print '\n'.join(mail.split("From")[-1].splitlines()[1:])
else:
    mail_num=int(sys.argv[1])

    print '\n'.join(['\n'.join(e.splitlines()[1:]) for e in \
                     mail.split("From")[-1:-2*mail_num:-2]][::-1])
