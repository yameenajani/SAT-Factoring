import requests
import json
import sys

if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    with open(input_filename) as fp:
        f = json.load(fp)
        link = f['link']
        res = requests.get(link).text.split("\n")[1:-2]
        res[0] = res[0][5:]
        res = "\n".join(res)
    with open(output_filename, "w") as file:
        file.write(res)
    
    with open(output_filename, "r") as file:
        lines = file.readlines()

    l = lines[9]
    l = l.split(" ")
    l[3] = str(int(l[3]) - 1)
    lines[9] = " ".join(l) + "\n"

    with open(output_filename, "w") as file:
        for line in lines:
            file.write(line)
