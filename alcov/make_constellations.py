import requests
import json
import time
import os
from tqdm import tqdm  # progress bar

# Global settings
BASE_URL = "https://lapis.cov-spectrum.org/open/v2"
HEADERS = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "application/json",
}
DATA_DIR = "data/constellations"
os.makedirs(DATA_DIR, exist_ok=True)

FAILED_FILE = "omitted_lineages.txt"
FAILED_LINEAGES = set()

# Load any previously failed lineages
if os.path.exists(FAILED_FILE):
    with open(FAILED_FILE, "r") as f:
        FAILED_LINEAGES = set(f.read().splitlines())
    print(f"🔄 Loaded {len(FAILED_LINEAGES)} previously failed lineages.")


def safe_get(url, lineage, retries=3, delay=5, timeout=20):
    """Make a GET request with retries and error handling. Logs failed lineages."""
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, headers=HEADERS, timeout=timeout)
            r.raise_for_status()
            data = r.json()

            # API-level error
            if "error" in data:
                detail = data["error"].get("detail", "Unknown error")
                print(f"❌ API error for {lineage}: {detail}")

                # If lineage is not valid, record & skip permanently
                if "not a valid lineage" in detail:
                    FAILED_LINEAGES.add(lineage)
                    return None

                # Otherwise retry if attempts remain
                elif attempt < retries:
                    time.sleep(delay)
                    continue
                else:
                    FAILED_LINEAGES.add(lineage)
                    return None

            return r
        except requests.exceptions.RequestException as e:
            print(f"⚠️ Attempt {attempt} for {lineage} at {url} failed: {e}")
            time.sleep(delay)

    FAILED_LINEAGES.add(lineage)
    return None


# Download master lineage list from cov-lineages pango-designation
input_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
response = requests.get(input_url)
if response:
    lines = response.text.splitlines()
    values = [line.split()[0] for line in lines[1:] if line.strip() and not line.startswith("*")]
    with open("lineages.txt", "w") as file:
        file.write("\n".join(values))
        file.close()
else:
    raise RuntimeError("Failed to download lineage list.")

# Read lineage list
with open("lineages.txt", "r") as file:
    lineages = file.read().splitlines()

# Filter out previously failed lineages
lineages = [l for l in lineages if l not in FAILED_LINEAGES]

# Iterate over lineages with progress bar
for lineage in tqdm(lineages, desc="Processing lineages", unit="lineage"):
    output_file = os.path.join(DATA_DIR, f"{lineage}.json")

    url1 = f"{BASE_URL}/sample/aggregated?dateFrom=2022-01-01&nextcladePangoLineage={lineage}"
    url3 = f"{BASE_URL}/sample/aggregated?dateFrom=2022-01-01&pangoLineage={lineage}"

    response1 = safe_get(url1, lineage)
    response3 = safe_get(url3, lineage)

    if not response1 or not response3:
        continue

    try:
        data1 = response1.json()
        data3 = response3.json()
        count1 = data1["data"][0]["count"] if data1["data"] else 0
        count3 = data3["data"][0]["count"] if data3["data"] else 0
    except Exception:
        FAILED_LINEAGES.add(lineage)
        continue

    nucleotide_mutations = []

    if count1 > 100:
        url2 = f"{BASE_URL}/sample/nucleotideMutations?nextcladePangoLineage={lineage}"
        response2 = safe_get(url2, lineage)
        if response2:
            data = response2.json()
            nucleotide_mutations = [
                m["mutation"] for m in data.get("data", [])
                if m.get("proportion", 0) > 0.90 and not m["mutation"].endswith("-")
            ]

    elif count3 > 100:
        url4 = f"{BASE_URL}/sample/nucleotideMutations?pangoLineage={lineage}"
        response4 = safe_get(url4, lineage)
        if response4:
            data = response4.json()
            nucleotide_mutations = [
                m["mutation"] for m in data.get("data", [])
                if m.get("proportion", 0) > 0.90 and not m["mutation"].endswith("-")
            ]

    if not nucleotide_mutations:
        continue

    json_data = {
        "label": lineage + "-like",
        "description": f"{lineage} lineage defining mutations",
        "sources": [],
        "tags": [lineage],
        "sites": nucleotide_mutations,
        "note": "Unique mutations for sublineage",
    }

    with open(output_file, "w") as file:
        json.dump(json_data, file)

# --- Write out all failed lineages ---
if FAILED_LINEAGES:
    with open(FAILED_FILE, "w") as f:
        f.write("\n".join(sorted(FAILED_LINEAGES)))
    print(f"\n❗ Saved {len(FAILED_LINEAGES)} failed lineages to {FAILED_FILE}")
else:
    print("\n✅ No failed lineages detected")
