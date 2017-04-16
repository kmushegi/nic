import os
import sys
import time

from twilio.rest import Client

accound_sid = "AC66228029c3056e9f58994d70688ce37b"
auth_token = "9a36b1a79ac0d0445b0003604a6f7047"

client = Client(accound_sid,auth_token)

num_files = len(os.listdir("stats"))
fifteenMinutesInSeconds = 15*60

while(num_files != 2187):
	num_files = len(os.listdir("stats"))
	message = client.messages.create(
			to="+12072086505",
			from_="+12077473393",
			body=("Currently " + str(num_files) + " files. That's " + str((float(num_files)/2187.0)*100)[:5] + "%"))
	print(message.sid)
	print("Now sleeping for 15 minutes")
	time.sleep(fifteenMinutesInSeconds)

if(num_files == 2187):
	message = client.messages.create(
			to="+12072086505",
			from_="+12077473393",
			body="Testing is done!")
	print(message.sid)

		