from IPython.display import display, HTML
import pandas as pd
import numpy, scipy
import matplotlib, matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import date
import seaborn as sns

pd.set_option('display.max_rows', 500)
version = %env WORKSPACE_CDR

condition_codes_icd9 = ['331.0', '331.11', '331.19 ', '331.82']

condition_codes_icd10 = ['G30.9', 'G31.01', 'G31.83']

drug_names = ['donepezil', 'aricept', 'galantamine', 'razadyne', 'reminyl', 'nivalin', 'rivastigmine', 'exelon', 
              'tacrine', 'cognex', 'memantine', 'namenda', 'axura', 'ebixa', 'akatinol']

condition_codes_icd9_subquery = ", ".join("'" + str(x) + "'" for x in condition_codes_icd9)
condition_codes_icd10_subquery = ", ".join("'" + str(x) + "'" for x in condition_codes_icd10)
condition_codes_icd9_subquery

condition_codes_icd9_subquery = ", ".join("'" + str(x) + "'" for x in condition_codes_icd9)
condition_codes_icd10_subquery = ", ".join("'" + str(x) + "'" for x in condition_codes_icd10)


Condition_Concept_IDs = pd.read_gbq("""
                                     
SELECT 
    concept_name as Condition_Name,
    concept_code as ICD_Code,
    concept_id as Concept_ID
FROM 
    `{0}.concept`
WHERE
    (vocabulary_id LIKE '%ICD9%' AND concept_code IN ({1}))
    OR (vocabulary_id LIKE '%ICD10%' AND concept_code IN ({2}))
GROUP BY
    Condition_Name,
    ICD_Code,
    Concept_ID
ORDER BY
    ICD_Code ASC

""".format(version, condition_codes_icd9_subquery, condition_codes_icd10_subquery), 
dialect = "standard")


Condition_Concept_IDs.shape

drug_list_subquery = " OR ".join("LOWER(c.concept_name) LIKE '%" + str(x) + "%'" for x in drug_names)

Drug_Concept_IDs = pd.read_gbq("""

SELECT
    DISTINCT c2.concept_name as Drug_Name,
    c2.concept_code as Concept_Code,
    c2.concept_id as Concept_ID
FROM
    `{0}.concept` c
    JOIN `{0}.concept_ancestor` ca
        ON c.concept_id = ca.ancestor_concept_id
    JOIN `{0}.concept` c2
        ON c2.concept_id = ca.descendant_concept_id
WHERE
    c.concept_class_id = 'Ingredient'
    AND ({1})

""".format(version, drug_list_subquery),
dialect = "standard")

Drug_Concept_IDs.shape

condition_concepts_subquery = ", ".join(str(x) for x in Condition_Concept_IDs['Concept_ID'])
condition_concepts_subquery
condition_concepts_subquery = ", ".join(str(x) for x in Condition_Concept_IDs['Concept_ID'])
#condition_concepts_subquery = ", ".join("'" + str(x) + "'" for x in Condition_Concept_IDs['ICD_Code'])

condition_PIDs = pd.read_gbq(f"""

SELECT
    person_id
FROM 
    `{version}.condition_occurrence`
WHERE 
    condition_source_concept_id IN ({condition_concepts_subquery})
GROUP BY
    person_id
HAVING
    COUNT(person_id) >= 5

""")['person_id'].tolist()

drug_concepts_subquery = ", ".join(str(x) for x in Drug_Concept_IDs['Concept_ID'])


drug_PIDs = pd.read_gbq(f"""

SELECT 
    DISTINCT p.person_id
FROM 
    `{version}.person` p
    LEFT JOIN `{version}.drug_exposure` d ON p.person_id = d.person_id
WHERE 
    drug_concept_id IN ({drug_concepts_subquery})

""")['person_id'].tolist()
cohort_PIDs = set(condition_PIDs).union(set(drug_PIDs))
print(
      "--------------\nFinal Results\n--------------\n\n"
    
      "Condition concepts - " + str(Condition_Concept_IDs.shape[0]) + 
      "\n" +
      "Drug concepts - " + str(Drug_Concept_IDs.shape[0]) + 
      "\n\n" +

      "Participants with condition - " + str(len(condition_PIDs)) + 
      "\n" +
      "Participants with drug exposure - " +  str(len(drug_PIDs)) + 
      "\n\n" +

      "Final cohort - " + str(len(cohort_PIDs))
     )



sns.set(style='darkgrid')
def round(arr):
    rounder = 20
    for x in range(len(arr)):
        if arr[x] < 20:
            arr[x] = (arr[x]//rounder + 1) * rounder
    return arr

cohort_subquery = ", ".join(str(x) for x in cohort_PIDs)


cohort = pd.read_gbq("""

SELECT
    DISTINCT person_id AS PERSON_ID,
    birth_datetime AS DATE_OF_BIRTH,
    c_race.concept_name AS RACE,
    c_sex.concept_name AS GENDER,
    c_ethn.concept_name AS ETHNICITY
FROM
    `{0}.person` p
    LEFT JOIN `{0}.concept` c_race
        ON p.race_concept_id = c_race.concept_id
    LEFT JOIN `{0}.concept` c_sex
        ON p.sex_at_birth_concept_id = c_sex.concept_id
    LEFT JOIN `{0}.concept` c_ethn
        ON p.ethnicity_concept_id = c_ethn.concept_id
WHERE person_id IN ({1})

""".format(version, cohort_subquery), dialect = "standard")


clean_cohort = cohort.rename(columns={'PERSON_ID':"Count",
                                      'DATE_OF_BIRTH':"Age",
                                      'RACE':"Race",
                                      'GENDER':"Sex at Birth",
                                      'ETHNICITY':"Hispanic"})


for row in (range(clean_cohort.shape[0])):
    for col in (range(clean_cohort.shape[1])):
        if clean_cohort.iloc[row,col] == "PMI: Skip":
            clean_cohort.iloc[row,col] = "Skip"
        if clean_cohort.iloc[row,col] in ["Not man only, not woman only, prefer not to answer, or skipped","Not male, not female, prefer not to answer, or skipped",
                                          "No matching concept",
                                          "None of these",
                                          "I prefer not to answer","What Race Ethnicity: Race Ethnicity None Of These","PMI: Prefer Not To Answer"]:
            clean_cohort.iloc[row,col] = "Unspecified"

            
for x in range(len(clean_cohort)):
    birth_date = clean_cohort.at[x,'Age']
    try:  
        birthday = birth_date.replace(year = date.today().year) 

    # raised when birth date is February 29 and the current year is not a leap year 
    except ValueError:  
        birthday = birth_date.replace(year = date.today().year, month = birth_date.month + 1, day = 1) 

    if birthday.date() > date.today(): 
        clean_cohort.at[x,'Age'] = date.today().year - birth_date.year - 1
    else: 
        clean_cohort.at[x,'Age'] = date.today().year - birth_date.year


bins = [18,26,36,46,56,66,76,86,100]
labels = ['18-25','26-35','36-45','46-55','56-65','66-75','76-85','86+']
clean_cohort['Age Group'] = pd.cut(clean_cohort['Age'], bins=bins, labels=labels, right=False)





