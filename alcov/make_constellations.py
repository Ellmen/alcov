# -*- coding: utf-8 -*-
"""
Created on Wed May 31 2023

@author: Jenn Knapp
email: jknapp@uwaterloo.ca
"""
"""
Purpose: Retrieve lineage-definining mutations from Cov-Spectrum for all lineages that have deposited sequences

Reguires:
lineages.txt (All lineages listed when searching CoV-Spectrum with the params World, AllTimes, and blank search bar)

"""
import requests
import json

# Download master lineage list from cov-lineages pango-designation (updated regularly as new lineages arise)
input_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
output_file = "lineages.txt"

response = requests.get(input_url)
lines = response.text.splitlines()

# Parse lineage list to exclude header, descriptions, and withdrawn lineages 
values = [line.split()[0] for line in lines[1:] if not line.startswith("*")]

# Write all lineage names to a file called lineages.txt
with open(output_file, "w") as file:
    file.write("\n".join(values))
    file.close()  # Close the file to ensure data is written to disk

base_url = "https://lapis.cov-spectrum.org/open/v1"

# Read lineage list
with open("lineages.txt", "r") as file:
    lineages = file.read().splitlines()

# Iterate over lineages
for lineage in lineages:
    output_file = "data/constellations/" + lineage + ".json"

# Construct URL
    url1 = base_url + "/sample/aggregated?nextcladePangoLineage=" + lineage

    try:
        # Send API request
        response = requests.get(url1)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the response as JSON
            data1 = response.json()

            # Extract the required data from the JSON response
            count = data1["data"][0]["count"]

            if count > 100:

                # Construct URL
                url2 = base_url + "/sample/nuc-mutations?nextcladePangoLineage=" + lineage

                try:
                    # Send API request
                    response = requests.get(url2)

                    # Check if the request was successful (status code 200)
                    if response.status_code == 200:
                        # Parse the response as JSON
                        data = response.json()

                        # Extract the required data from the JSON response
                        nucleotide_mutations = [
                            mutation["mutation"] for mutation in data["data"]
                            if mutation["proportion"] > 0.90 and not mutation["mutation"].endswith("-")
                        ]

                        # Create the final JSON structure
                        json_data = {
                            "label": lineage + "-like",
                            "description": f"{lineage} lineage defining mutations",
                            "sources": [],
                            "tags": [lineage],
                            "sites": nucleotide_mutations,
                            "note": "Unique mutations for sublineage",
                            "rules": {
                                "default": {
                                    "min_alt": "",
                                    "max_ref": ""
                                },
                                "Probable": {
                                    "min_alt": "",
                                    "max_ref": ""
                                }
                            }
                        }

                        # Write data to JSON file
                        with open(output_file, "w") as file:
                            json.dump(json_data, file)
                    else:
                        print("Error occurred while fetching data for lineage:", lineage)
         
                except requests.RequestException as e:
                        print("An error occurred during the API request:", str(e))

    except requests.RequestException as e:
        print("An error occurred during the API request:", str(e))

