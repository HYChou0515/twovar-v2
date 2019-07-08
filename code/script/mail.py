# 1. Save the encrypted_hash from encrypt.py to 'gmail.pub'
# 2. Change 'samepl@gmail.com' to your email address
# 3. Usage: see send_mail.py as an example

import smtplib
from encrypt import check_passwd
from getpass import getpass
from os.path import abspath, dirname
import os,sys
import logging
import traceback


def _get_logger(parent_logger, name):
	if parent_logger is None:
		plogger = logging.getLogger('mail.'+name)
		plogger.setLevel(logging.WARN)
		ch = logging.StreamHandler()
		ch.setLevel(logging.WARN)
		formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
		ch.setFormatter(formatter)
		plogger.addHandler(ch)
	else:
		plogger = parent_logger.getChild(name)
	return plogger

class Mail(object):
	def __init__(self, mailaddr='sample@gmail.com', mail_pub_f='gmail.pub', logger=None):
		if logger is None:
			logger = self.__get_default_logger__()
		# get pub key from mail_pub_f
		with open(os.path.join(dirname(abspath(__file__)), mail_pub_f)) as f:
			mail_pub = f.read()
		# get passwd
		if sys.stdin.isatty():
			self.mail_passwd = getpass('gmail password: ')
		else:
			self.mail_passwd = sys.stdin.readline().rstrip()
		# check if passwd match pub key
		if not check_passwd(self.mail_passwd, mail_pub):
			raise ValueError('Wrong Password')
		self.logger = _get_logger(logger, 'mail')
		self.mailaddr = mailaddr
		self.msgs = []
	def __get_default_logger__(self):
		default_logger = logging.getLogger()
		default_logger.setLevel(logging.WARN)

		handler = logging.StreamHandler(sys.stdout)
		handler.setLevel(logging.WARN)
		formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
		handler.setFormatter(formatter)
		default_logger.addHandler(handler)

		return default_logger

	def append_msg(self, msg):
		self.logger.info('append message'+msg)
		self.msgs.append(msg)
	def send(self, title):
		try:
			server = smtplib.SMTP('smtp.gmail.com', 587)
			server.ehlo()
			server.starttls()
			server.login(self.mailaddr, self.mail_passwd)
		except:
			self.logger.error(traceback.format_exc())
		msg = "\r\n".join([
			"From: "+self.mailaddr,
			"To: "+self.mailaddr,
			"Subject: "+title,
			"",
			'\n'.join(self.msgs),
		])
		try:
			server.sendmail(self.mailaddr, self.mailaddr, msg)
			server.quit()
		except:
			self.logger.error(traceback.format_exc())
