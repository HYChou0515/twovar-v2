#!/usr/bin/env python
# Usage: ./send_mail.py mail_address pubkey_file
#   and enter your passwd
# Other Usage: echo "password" | ./send_mail.py mail_address pubkey_file

from mail import Mail

def main():
	import sys
	def exit_with_help():
		print('Usage: %s mail_address pubkey_file')
		sys.exit(1)
	argv = sys.argv
	if len(argv) != 3:
		exit_with_help()
	mailaddr = argv[1]
	pub_f = argv[2]

	mail = Mail(mailaddr, pub_f)
	mail.append_msg('This is containt')
	mail.send('This is mail title')

if __name__ == '__main__':
	main()
