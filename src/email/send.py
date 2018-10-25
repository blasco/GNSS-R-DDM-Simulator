#!/usr/bin/env python

def send_simple_message():
    return requests.post(
        "https://api.mailgun.net/v3/sandbox0ef7a74039aa4f4e8e7afd631f0e0dda.mailgun.org/messages",
        auth=("api", ""),
        data={"from": "Ecobloom <juan@ecobloom.se>",
              "to": "Juan Blasco <blascoburguillos@gmail.com>",
              "subject": "Hello Juan Blasco",
              "text": "Congratulations Juan Blasco, you just sent an email with Mailgun!  You are truly awesome!"})
