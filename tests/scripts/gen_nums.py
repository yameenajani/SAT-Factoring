#!/usr/bin/python3

from sage.all import ZZ, random_prime
import sys
import json
import os.path
from bs4 import BeautifulSoup
import urllib.request

def generate_test_case(bitsize):
    p = ZZ(0)
    q = ZZ(0)
    while (p-1)%3 != 1 or (q-1)%3 != 1:
        p = random_prime(2**bitsize-1,False,2**(bitsize-1))
        q = random_prime(2**bitsize-1,False,2**(bitsize-1))
    n = p*q
    return p, q, n

def main():
    bitsizes = [int(sys.argv[1])]
    num_instances = int(sys.argv[2])
    for bitsize in bitsizes:
        for i in range(1, num_instances+1):
            # Uncomment below code to update links to the insatance
            ''' if os.path.exists("../data/{}/data_{}.json".format(bitsize, i)):
                with open("../data/{}/data_{}.json".format(bitsize, i), "r") as f:
                    data = json.load(f)
                    p = data['p']
                    q = data['q']
                    n = data['n']
            else: '''
            p, q, n = generate_test_case(bitsize)
            weblink = "https://cgi.luddy.indiana.edu/~sabry/cnf.cgi?factor={}&Adder=nbit&Multiplier=recursive".format(n)
            html_page = urllib.request.urlopen(weblink)
            soup = BeautifulSoup(html_page, "html.parser")
            for link in soup.findAll('a'):
                l = link.get('href')
            data = {
                "p": int(p) if p<q else int(q),
                "q": int(q) if p<q else int(p),
                "n": int(n),
                "link": "https://cgi.luddy.indiana.edu/~sabry/" + "{}".format(l)
            }

            json_obj = json.dumps(data, indent=4)

            if not (os.path.exists("../data/{}".format(bitsize))):
                os.makedirs("../data/{}".format(bitsize))
            # Even if the file already exists, the weblink may need to be updated
            with open("../data/{}/data_{}.json".format(bitsize, i), "w") as f:
                f.write(json_obj)

if __name__ == '__main__':
    main()
