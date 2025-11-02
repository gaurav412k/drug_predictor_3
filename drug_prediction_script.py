import subprocess
import sys
import os
import re
import json
from pathlib import Path

# List of required packages
required_packages = [
    'pandas',
    'sklearn'
]

# Function to install missing packages
def install_package(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Check and install required packages
for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        if package == "sklearn":
            print(f"Package '{package}' not found. Installing 'scikit-learn'...")
            install_package("scikit-learn==1.5.2")
        else:
            print(f"Package '{package}' not found. Installing...")
            install_package(package)

# Now proceed with the main script
import pickle
import pandas as pd # type: ignore

with open(r'model_data.json', 'r', encoding='utf-8') as f:
    json_data = json.load(f)


# Load the saved models and encoders from the pickle file
with open("models_and_encoders.pkl", 'rb') as f:
    models_and_encoders = pickle.load(f)

# Extract models and encoders from the dictionary
model_1 = models_and_encoders['model_1']
model_2 = models_and_encoders['model_2']
target_1_encoder = models_and_encoders['target_1_encoder']
target_2_encoder = models_and_encoders['target_2_encoder']
# Converting the data the from the input file for the model,like lowercasing the bacteria names and converting the gene name to resistance names
# This is also a checkpoint to check if the input data is correct or not
def convertingFeatures(df,json_data=json_data):
    feature_type = df["type"].strip().lower()
    feature_input = df["input_feature"]

    if feature_type == "organism":
        organism = feature_input.strip().lower()
        if organism in json_data["all_features"]:
            return organism
        else:
            print(f"[Warning] '{organism}' is not a recognized organism. Please check spelling.")
            return None
    elif feature_type == "genes":
        gene_list = [g.strip().lower() for g in feature_input.split(",")]
        antibiotics = [json_data["gene_to_antibiotic"][g] for g in gene_list if g in json_data["gene_to_antibiotic"]]

        if antibiotics:
            return antibiotics[0]  # Return the first mapped antibiotic
        else:
            print(f"[Warning] No valid genes found in: {feature_input}. Please verify gene names.")
            return None
    else:
        print(f"[Warning] Unknown feature type: {feature_type}")
        return None

#This function gets the proper name from the predicted data
def get_drugs_name(text,json_data=json_data['synonyms']):
    reverse_lookup = {}
    for canon, syn_list in json_data.items():
        for s in syn_list:
            reverse_lookup[s.lower()] = canon

    escaped_syns = [re.escape(s) for s in reverse_lookup.keys()]
    pattern = re.compile(r'\b(' + '|'.join(escaped_syns) + r')\b', flags=re.IGNORECASE)

    drugs_found = { reverse_lookup[match.group(1).lower()] 
            for match in pattern.finditer(text) }

    return list(drugs_found)

#this function gets the score for each drug
def get_adverse_reaction_scores(drug_input,json_data=json_data["drug_score"]):
    drug_input = drug_input.split(",")  # Split input if there are multiple drugs
    for drug in drug_input:
        drug = drug.strip()  # Remove any leading/trailing spaces
        if drug in json_data:
            score_str = json_data[drug]
            # Handle case when the score is 'not applicable' or '-'
            if score_str == '-':
                return 1  # Return 1 for "not applicable"
            score_list = score_str.split(',')
            try:
                # Convert the first score to an integer
                score = int(score_list[0])
                return score
            except ValueError:
                # In case there's an issue converting to int, return default score
                return 1
    return 1  # If no score found, return default score


#This function is use to get the resistance score from the given json file
def get_antimicrobial_resistance_score(resistances,json_data=json_data["resistance_score"]):
    total_score = sum(json_data.get(res.lower(), 0) for res in resistances)
    return total_score

#This fucntion is to get the side effects and interactions from the json data for the primary drug
def get_side_effects_interactions(primary_drug,data,json_data=json_data["drug_effects"]):
    return json_data[primary_drug][data]

# This function adds the string of simple or complicated uti based on the dosage provided in the primary drug
def detect_uti_type(text):
    text_lower = text.lower()

    if 'doses' in text_lower:
        return " for possible complicated UTI *"
    elif 'dose' in text_lower:
        return " for possible simple UTI *"


    matches = re.findall(r'(\d+)(?:-(\d+))?\s*days?', text_lower)
    for match in matches:
        start = int(match[0])
        end = int(match[1]) if match[1] else start

        if match[1]:  # range like 5-7
            if start < 7:
                return " for possible simple UTI *"
            else:
                return " for possible complicated UTI *"
        else:  # single number
            if start >= 7:
                return " for possible complicated UTI *"
            else:
                return " for possible simple UTI *"

    return ""  # Default fallback

def replace_drug_synonyms_with_canonical_names(text, json_data=json_data['synonyms']):
    # Build reverse lookup from synonym to canonical name
    reverse_lookup = {}
    for canon, syn_list in json_data.items():
        for s in syn_list:
            reverse_lookup[s.lower()] = canon

    # Build regex pattern
    escaped_syns = sorted([re.escape(s) for s in reverse_lookup.keys()], key=len, reverse=True)
    pattern = re.compile(r'\b(' + '|'.join(escaped_syns) + r')\b', flags=re.IGNORECASE)

    # Replacement function
    def replacer(match):
        matched_text = match.group(0)
        canonical = reverse_lookup.get(matched_text.lower())
        return canonical if canonical else matched_text

    # Replace and return modified text
    return pattern.sub(replacer, text)

# this function is to prepare the predicted data in the final form
def preparing_finalData(primary_drugs, secondary_drugs, data_df):
    annotation = detect_uti_type(primary_drugs)
    output_dict = {}

    # Clean and extract drug names
    clean_primary_drug = primary_drugs.strip("[]'").replace("| ", "").strip()
    primary_names = get_drugs_name(primary_drugs)
    secondary_names = get_drugs_name(secondary_drugs)

    # Populate output
    output_dict["Primary Drug"] = f"{replace_drug_synonyms_with_canonical_names(clean_primary_drug).replace("| ", "").strip()}{annotation}"
    output_dict["Primary Drug Score"] = [f"{drug} : {get_adverse_reaction_scores(drug)}" for drug in primary_names]
    output_dict["Primary Drug Side Effects"] = [f"{drug} : {get_side_effects_interactions(drug, 'side effects')}" for drug in primary_names]
    output_dict["Primary Drug Interactions"] = [f"{drug} : {get_side_effects_interactions(drug, 'interactions')}" for drug in primary_names]
    output_dict["Primary Drug Dose Adj"] = [f"{drug} : {get_side_effects_interactions(drug, 'dose adj')}" for drug in primary_names]

    output_dict["Secondary Drug"] = replace_drug_synonyms_with_canonical_names(secondary_drugs)
    output_dict["Secondary Drug Score"] = [f"{drug} : {get_adverse_reaction_scores(drug)}" for drug in secondary_names]

    resistances = data_df[data_df["type"] == "genes"]["finalFeature"].str.title().tolist()
    output_dict["Resistance"] = resistances
    output_dict["Resistance Score"] = get_antimicrobial_resistance_score(resistances)

    return output_dict


data_df = pd.read_csv(sys.argv[1],dtype={"ct_mean" : float})
data_df = data_df.dropna(subset=["type","input_feature","ct_mean"])
data_df["finalFeature"] = data_df.apply(convertingFeatures,axis=1)
# data_df
input_features = dict(
    zip(
        data_df.loc[data_df["finalFeature"].notnull(), "finalFeature"],
        data_df.loc[data_df["finalFeature"].notnull(), "ct_mean"]
    )
)
full_input = {feat: input_features.get(feat, 0.0) for feat in json_data["all_features"]}
print("==========================================================================================================")
print(f"Input Values are:\n{"\n".join(f"{key}: {val}" for key, val in input_features.items())}")
print("==========================================================================================================")
input_df = pd.DataFrame([full_input])
# 🔹 Predict target_1
y1_pred = model_1.predict(input_df)
input_df['target_1_encoded'] = y1_pred

# 🔹 Predict target_2
y2_pred = model_2.predict(input_df)

# 🔹 Decode predictions
decoded_y1 = target_1_encoder.inverse_transform(y1_pred)[0]
decoded_y2 = target_2_encoder.inverse_transform(y2_pred)[0]
output_dict = preparing_finalData(decoded_y1,decoded_y2,data_df)
output_dict
print("\n🎯 Predicted Results:")
print(f"🔹 Primary Drug: {output_dict['Primary Drug']}")
print(f"🔹 Primary Drug Score: {output_dict['Primary Drug Score']}")

print(f"🔹 Primary Drug Side Effects: {output_dict['Primary Drug Side Effects']}")
print(f"🔹 Primary Drug Interactions: {output_dict['Primary Drug Interactions']}")
print(f"🔹 Primary Drug Dose Adj: {output_dict['Primary Drug Dose Adj']}")

print(f"🔹 Secondary Drugs:\n{output_dict['Secondary Drug']}")
print(f"🔹 Secondary Drugs Score: {output_dict['Secondary Drug Score']}")

print(f"🔹 Resistance: {output_dict['Resistance']}")
print(f"🔹 Resistance Score: {output_dict['Resistance Score']}")

finalDf = pd.DataFrame([output_dict]).transpose().reset_index()
finalDf.rename(columns={0 : "Values","index" : "Category"},inplace=True)
# finalDf

finalDf.to_excel(r"Predicted_Drug_Output.xlsx",index=False)
print(f"\n\nThe output is saved with Drug_Output.xlsx name in {os.path.dirname(os.path.realpath(__file__))}")