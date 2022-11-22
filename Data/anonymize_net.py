import pandas
import json
import pandas
import os


file_path = os.path.join()

with open (file_path) as f:
    data=json.load(f)
#read in network#
data=data["net_name"]['interaction']
json_object = json.dumps(data, indent=4)
with open ("Data/disease_net.json", "w") as f:
    f.write(json_object)
    f.close()

