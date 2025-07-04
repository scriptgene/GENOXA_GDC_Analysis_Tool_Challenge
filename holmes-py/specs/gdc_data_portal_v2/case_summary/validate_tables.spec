# Case Summary - Validate Tables
Date Created        : 09/05/2024
Version			        : 1.0
Owner		            : GDC QA
Description		      : Validate All Tables and Sections in Case Summary Page
Test-Case           : PEAR-2139, PEAR-464

tags: gdc-data-portal-v2, case-summary, clinical-biospecimen-download, regression

## Navigate to Case Summary Page: TCGA-13-0920
* On GDC Data Portal V2 app
* Quick search for "85a85a11-7200-4e96-97af-6ba26d680d59" and go to its page
## Summary Table
* Verify the table "Summary Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Case UUID                              |
    |Case ID                                |
    |Project Name                           |
    |Disease Type                           |
    |Program                                |
    |Primary Site                           |
    |Images                                 |
    |85a85a11-7200-4e96-97af-6ba26d680d59   |
    |TCGA-OV                                |
    |Cystic, Mucinous and Serous Neoplasms  |
    |Ovary                                  |
* Select "Add Remove Files" on Case Summary page
* Is modal with text "Added 2 files to the cart." present on the page and "Remove Modal"
* Select "Add Remove Files" on Case Summary page
* Is modal with text "Removed 2 files from the cart." present on the page and "Remove Modal"
* Select button "View Slide Images"
* Select "Details" on the Image Viewer page
* Verify details fields and values
  |field_name                       |value                                  |
  |---------------------------------|---------------------------------------|
  |File_id                          |010fc715-0098-4967-a36f-a8922c03e40e   |
  |Submitter_id                     |TCGA-13-0920-01A-01-BS1                |
  |Slide_id                         |b7b9daea-d790-433c-93b2-db5e26827f30   |

## Navigate Back To Case Summary Page: TCGA-13-0920
* Quick search for "85a85a11-7200-4e96-97af-6ba26d680d59" and go to its page
## Data Category Table
* Verify the table "Data Category Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Data Category                          |
    |Biospecimen                            |
    |Clinical                               |
    |Copy Number Variation                  |
    |DNA Methylation                        |
    |Proteome Profiling                     |
    |Sequencing Reads                       |
    |Simple Nucleotide Variation            |
    |Structural Variation                   |
    |Transcriptome Profiling                |

## Experimental Strategy Table
* Verify the table "Experimental Strategy Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Experimental Strategy                  |
    |Expression Array                       |
    |Genotyping Array                       |
    |Methylation Array                      |
    |miRNA-Seq                              |
    |Reverse Phase Protein Array            |
    |Tissue Slide                           |
    |WGS                                    |
    |WXS                                    |

## Navigate Back To Case Summary Page: TARGET-20-PASGWS
* Quick search for "a51ed206-81e4-42d8-a50d-d19cd66f5dca" and go to its page
## Clinical Table - No Data
* Select tab "Diagnoses Treatments" on Case Summary page
* Is text "No Diagnoses Found." present on the page

* Select tab "Family Histories" on Case Summary page
* Is text "No Family Histories Found." present on the page

* Select tab "Exposures" on Case Summary page
* Is text "No Exposures Found." present on the page

* Select tab "Followups Molecular Tests" on Case Summary page
* Is text "No Follow Ups Found." present on the page

* Select tab "Demographic" on Case Summary page
* Is text "No Demographic Found." present on the page

## Navigate to Case Summary Page: C3N-03018
* Quick search for "5fd8d17e-57d2-4270-8728-de259ff6b2fe" and go to its page
## Clinical Table - No Molecular Tests, No Other Clinical Attributes, No Treatments
* Select tab "Followups Molecular Tests" on Case Summary page
* Is text "No Molecular Tests Found." present on the page
* Is text "No Other Clinical Attributes Found." present on the page

* Select tab "Diagnoses Treatments" on Case Summary page
* Is text "No Treatments Found." present on the page

## Clinical Table - Exposure Tab
* Select tab "Exposures" on Case Summary page
* Verify the table "Clinical Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Exposure ID                            |
    |Exposure UUID                          |
    |Alcohol History                        |
    |Alcohol Intensity                      |
    |Tobacco Smoking Status                 |
    |Pack Years Smoked                      |
    |C3N-03018-EXP                          |
    |f1416551-58ec-43bb-a1f0-78658dae0794   |
    |Yes                                    |
    |Occasional Drinker                     |
    |Current Reformed Smoker for > 15 yrs   |
## Clinical Table - TSV Download: Exposure, Pathological Detail
* Download "TSV" from "Case Summary Clinical Table"
* Read file content from compressed "TSV from Case Summary Clinical Table"
* Verify that "TSV from Case Summary Clinical Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |5fd8d17e-57d2-4270-8728-de259ff6b2fe   |
    |exposures.cigarettes_per_day           |
    |exposures.exposure_type                |
    |Tobacco                                |
    |exposures.respirable_crystalline_silica_exposure |
    |exposures.tobacco_smoking_status       |
    |Current Reformed Smoker for > 15 yrs   |
    |diagnosis_id                           |
    |26a5ccaa-d86d-4283-9a00-39ed229586f9   |
    |extracapsular_extension_present        |
    |pathology_details.margin_status        |
    |pathology_details.submitter_id         |
    |3.5                                    |
* Verify that "TSV from Case Summary Clinical Table" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |tumor_stage                            |
    |marijuana_use_per_week                 |
    |smokeless_tobacco_quit_age             |
    |tobacco_use_per_day                    |

## Clinical Table - Download JSON: Exposure, Pathological Detail
* Download "JSON" from "Case Summary Clinical Table"
* Read from "JSON from Case Summary Clinical Table"
* Verify that "JSON from Case Summary Clinical Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |Tobacco                                |
    |Current Reformed Smoker for > 15 yrs   |
    |1968                                   |
    |Stage I                                |
    |Renal cell carcinoma, NOS              |
    |30196a26-d73e-428c-9f7a-707e6a2c7241   |
    |C3N-03018-PATH                         |
    |Kidney, NOS                            |

## Navigate to Case Summary Page: MMRF_2081
* Quick search for "3d5ebf3f-0cbd-458a-820d-65652e9682d7" and go to its page
## Clinical Table - Demographic
* Select tab "Demographic" on Case Summary page
* Verify the table "Clinical Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Demographic ID                         |
    |Demographic UUID                       |
    |Ethnicity                              |
    |Gender                                 |
    |Race                                   |
    |Days To Birth                          |
    |Days To Death                          |
    |Vital Status                           |
    |MMRF_2081_demographic1                 |
    |857aa4ee-8758-4297-9f67-539ff420c1cd   |
    |not hispanic or latino                 |
    |female                                 |
    |white                                  |
    |-69 years 151 days                     |
    |--                                     |
    |Alive                                  |

## Clinical Table - Diagnoses and Treatments
* Select tab "Diagnoses Treatments" on Case Summary page
* Verify the table "Clinical Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Diagnosis ID                           |
    |Diagnosis UUID                         |
    |Classification Of Tumor                |
    |Age At Diagnosis                       |
    |Days To Last Follow Up                 |
    |Days To Last Known Disease Status      |
    |Days To Recurrence                     |
    |Morphology                             |
    |Primary Diagnosis                      |
    |Prior Malignancy                       |
    |Synchronous Malignancy                 |
    |Progression Or Recurrence              |
    |Site Of Resection Or Biopsy            |
    |Tissue Or Organ Of Origin              |
    |Tumor Grade                            |
    |MMRF_2081_diagnosis1                   |
    |0ad0643b-8364-4ea3-a600-d03c0be2e6fd   |
    |69 years 151 days                      |
    |Unknown tumor status                   |
    |9732/3                                 |
    |Multiple myeloma                       |
    |Therapeutic Agents                     |
    |Treatment Intent Type                  |
    |Treatment or Therapy                   |
    |Days to Treatment Start                |
    |MMRF_2081_treatment2                   |
    |8af43e39-64d4-47df-b555-53c899377182   |
    |Dexamethasone                          |
    |180 days                               |

## Clinical Table - Family Histories
* Select tab "Family Histories" on Case Summary page
* Verify the table "Clinical Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Family History ID                      |
    |Family History UUID                    |
    |Relationship Age At Diagnosis          |
    |Relationship Gender                    |
    |Relationship Primary Diagnosis         |
    |Relationship Type                      |
    |Relative With Cancer History           |
    |MMRF_2081_family_history2              |
    |09d90ab1-aeb5-4bd8-9f1a-e2aa499251f0   |
    |male                                   |
    |Unknown                                |
    |Father                                 |

## Clinical Table - Followups, Molecular Tests, Other Clinical Attributes
* Select tab "Followups Molecular Tests" on Case Summary page
* Verify the table "Clinical Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Follow Up ID                           |
    |Follow Up UUID                         |
    |Days To Follow Up                      |
    |Progression Or Recurrence Type         |
    |Disease Response                       |
    |ECOG Performance Status                |
    |Karnofsky Performance Status           |
    |Progression Or Recurrence Anatomic Site|
    |MMRF_2081_followup1                    |
    |b2bcc701-7758-44fa-a118-6730b14dd542   |
    |-4 days                                |
    |0                                      |
    |Gene Symbol                            |
    |Second Gene Symbol                     |
    |Molecular Analysis Method              |
    |Laboratory Test                        |
    |Test Value                             |
    |Test Result                            |
    |Test Units                             |
    |Biospecimen Type                       |
    |Variant Type                           |
    |Chromosome                             |
    |AA Change                              |
    |Antigen                                |
    |Mismatch Repair Mutation               |
    |MMRF_2081_molecular_test17             |
    |a432493e-dd05-41db-bf9b-5c0ac46c4429   |
    |Absolute Neutrophil                    |
    |2.8                                    |
    |ukat/L                                 |
    |Timepoint Category                     |
    |Nononcologic Therapeutic Agents        |
    |Treatment Frequency                    |
    |Weight                                 |
    |Height                                 |
    |BMI                                    |
    |58                                     |
    |150                                    |
    |a5de9cbc-4033-4d39-b543-e3ad87ef0d78   |
    |MMRF_2081_oca1                         |

## Clinical Table - TSV Download: Follow Up, Family History, Clinical
* Download "TSV" from "Case Summary Clinical Table"
* Read file content from compressed "TSV from Case Summary Clinical Table"
* Verify that "TSV from Case Summary Clinical Table" has expected information
    |required_info                                    |
    |-------------------------------------------------|
    |follow_ups.cause_of_response                     |
    |follow_ups.barretts_esophagus_goblet_cells_present|
    |follow_ups.hysterectomy_margins_involved         |
    |follow_ups.recist_targeted_regions_number        |
    |Lactate Dehydrogenase                            |
    |molecular_tests.molecular_analysis_method        |
    |pathology_details.bone_marrow_malignant_cells    |
    |pathology_details.rhabdoid_present               |
    |pathology_details.zone_of_origin_prostate        |
    |family_histories.relative_smoker                 |
    |family_histories.relatives_with_cancer_history_count|
    |family_histories.relationship_gender             |
    |3d5ebf3f-0cbd-458a-820d-65652e9682d7             |
    |MMRF_2081                                        |
    |pathology_details.additional_pathology_findings  |
    |treatments.treatment_type                        |
    |MMRF-COMMPASS                                    |
    |other_clinical_attributes.aids_risk_factors      |
    |other_clinical_attributes.days_to_risk_factor    |
    |other_clinical_attributes.hysterectomy_margins_involved|
    |other_clinical_attributes.timepoint_category     |
    |other_clinical_attributes.undescended_testis_corrected_laterality|
    |other_clinical_attributes.weight                 |
    |MMRF_2081_oca1                                   |
    |a5de9cbc-4033-4d39-b543-e3ad87ef0d78             |
    |cases.case_id                                    |
    |cases.consent_type                               |
    |cases.days_to_consent                            |
    |cases.days_to_lost_to_followup                   |
    |cases.disease_type                               |
    |cases.index_date                                 |
    |cases.lost_to_followup                           |
    |cases.primary_site                               |
    |cases.submitter_id                               |


* Verify that "TSV from Case Summary Clinical Table" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |tumor_stage                            |
    |marijuana_use_per_week                 |
    |smokeless_tobacco_quit_age             |
    |tobacco_use_per_day                    |
    |diagnoses.anaplasia_present            |
    |diagnoses.anaplasia_present_type       |
    |diagnoses.breslow_thickness            |
    |diagnoses.circumferential_resection_margin|
    |diagnoses.greatest_tumor_dimension     |
    |diagnoses.gross_tumor_weight           |
    |diagnoses.largest_extrapelvic_peritoneal_focus|
    |diagnoses.lymph_node_involved_site     |
    |diagnoses.lymph_nodes_positive         |
    |diagnoses.lymph_nodes_tested           |
    |diagnoses.lymphatic_invasion_present   |
    |diagnoses.non_nodal_regional_disease   |
    |diagnoses.non_nodal_tumor_deposits     |
    |diagnoses.percent_tumor_invasion       |
    |diagnoses.perineural_invasion_present  |
    |diagnoses.peripancreatic_lymph_nodes_positive|
    |diagnoses.peripancreatic_lymph_nodes_tested|
    |diagnoses.transglottic_extension       |
    |diagnoses.tumor_largest_dimension_diameter|
    |diagnoses.tumor_stage                  |
    |diagnoses.vascular_invasion_present    |
    |diagnoses.vascular_invasion_type       |
    |exposures.bmi                          |
    |exposures.height                       |
    |exposures.marijuana_use_per_week       |
    |exposures.smokeless_tobacco_quit_age   |
    |exposures.tobacco_use_per_day          |
    |exposures.weight                       |

## Clinical Table - Download JSON: Follow Up, Family History, Clinical
* Download "JSON" from "Case Summary Clinical Table"
* Read from "JSON from Case Summary Clinical Table"
* Verify that "JSON from Case Summary Clinical Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |3d5ebf3f-0cbd-458a-820d-65652e9682d7   |
    |MMRF_2081                              |
    |MMRF_2081_family_history2              |
    |MMRF_2081_followup7                    |
    |MMRF-COMMPASS                          |
    |0ad0643b-8364-4ea3-a600-d03c0be2e6fd   |
    |not hispanic or latino                 |
    |other_clinical_attributes              |
    |other_clinical_attribute_id	          |
    |a5de9cbc-4033-4d39-b543-e3ad87ef0d78   |
    |MMRF_2081_oca1                         |
    |58                                     |
    |2024-07-08T14:50:58.508358-05:00       |
    |150                                    |
    |height                                 |
    |weight                                 |
    |Plasma Cell Tumors                     |


## Clinical Table - Validate JSON File Fields
  |field_name                               |
  |-----------------------------------------|
  |disease_type                             |
  |index_date                               |
  |primary_site                             |
  |case_id	                                |
  |submitter_id	                            |
  |family_histories.relative_with_cancer_history|
  |project.project_id                       |
* Verify that the "JSON from Case Summary Clinical Table" has <field_name> for each object


## Navigate to Case Summary Page: TCGA-GV-A3QI
* Quick search for "06f936e9-5a90-40d3-b91a-713f2b4e6e11" and go to its page
## Biospecimen Table - Expand/Collapse Button
* Select "Expand All"
* Verify the table "Biospecimen Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |samples                                |
    |TCGA-GV-A3QI-10A                       |
    |portions                               |
    |TCGA-GV-A3QI-10A-01                    |
    |analytes                               |
    |TCGA-GV-A3QI-10A-01                    |
    |aliquots                               |
    |TCGA-GV-A3QI-10A-01D-A21Y-01           |
    |TCGA-GV-A3QI-10A-01D-A21Z-08           |
    |slides                                 |
    |TCGA-GV-A3QI-01A-01-TSA                |
* Select "Collapse All"
* Pause "2" seconds
* Verify the table "Biospecimen Case Summary" is not displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |TCGA-GV-A3QI-10A-01                    |
    |analytes                               |
    |TCGA-GV-A3QI-10A-01                    |
    |aliquots                               |
    |TCGA-GV-A3QI-10A-01D-A21Y-01           |
    |TCGA-GV-A3QI-10A-01D-A21Z-08           |
    |slides                                 |
    |TCGA-GV-A3QI-01A-01-TSA                |
    |TCGA-GV-A3QI-01A-11D-A21Z-08           |
    |TCGA-GV-A3QI-01Z-00-DX1                |

## Biospecimen Table - Samples
* Enter "074df1cb-b6b7-4d3f-9c0a-fe523b1854e0" in the text box "Biospecimen Search Bar"
* Verify the table "Selection Information Biospecimen" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Sample ID                              |
    |Sample UUID                            |
    |Tissue Type                            |
    |Tumor Descriptor                       |
    |Specimen Type                          |
    |Preservation Method                    |
    |Tumor Code ID                          |
    |Shortest Dimension                     |
    |Intermediate Dimension                 |
    |Longest Dimension                      |
    |Pathology Report UUID                  |
    |Current Weight                         |
    |Initial Weight                         |
    |Time Between Clamping And Freezing     |
    |Time Between Excision And Freezing     |
    |Days To Sample Procurement             |
    |Freezing Method                        |
    |Days To Collection                     |
    |Portions                               |
    |TCGA-GV-A3QI-10A                       |
    |074df1cb-b6b7-4d3f-9c0a-fe523b1854e0   |
    |Peripheral Blood NOS                   |
    |Normal                                 |
    |Not Applicable                         |
    |105 days                               |

## Biospecimen Table - Portions
* Enter "1d4435fe-3a24-426d-9b49-d800551fc4a0" in the text box "Biospecimen Search Bar"
* Verify the table "Selection Information Biospecimen" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Portion ID                             |
    |Portion UUID                           |
    |Portion Number                         |
    |Weight                                 |
    |Is Ffpe	                              |
    |Analytes                               |
    |Slides                                 |
    |TCGA-GV-A3QI-01A-11                    |
    |1d4435fe-3a24-426d-9b49-d800551fc4a0   |
    |false                                  |
    |--                                     |

## Biospecimen Table - Slides
* Enter "88217dc9-06a0-4839-8269-85ce0798ef89" in the text box "Biospecimen Search Bar"
* Verify the table "Selection Information Biospecimen" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Slide ID                               |
    |Slide UUID                             |
    |Percent Tumor Nuclei                   |
    |Percent Monocyte Infiltration          |
    |Percent Normal Cells                   |
    |Percent Stromal Cells                  |
    |Percent Eosinophil Infiltration        |
    |Percent Lymphocyte Infiltration        |
    |Percent Neutrophil Infiltration        |
    |Section Location                       |
    |Percent Granulocyte Infiltration       |
    |Percent Necrosis                       |
    |Percent Inflam Infiltration            |
    |Number Proliferating Cells             |
    |Percent Tumor Cells                    |
    |Slide Image                            |
    |TCGA-GV-A3QI-01A-01-TSA                |
    |88217dc9-06a0-4839-8269-85ce0798ef89   |
    |TOP                                    |
    |100                                    |
* Select button "Add Remove Cart Biospecimen"
* Is modal with text "Added TCGA-GV-A3QI-01A-01-TSA.88217DC9-06A0-4839-8269-85CE0798EF89.svs to the cart." present on the page and "Remove Modal"
* Select button "Add Remove Cart Biospecimen"
* Is modal with text "Removed TCGA-GV-A3QI-01A-01-TSA.88217DC9-06A0-4839-8269-85CE0798EF89.svs from the cart." present on the page and "Remove Modal"
* Select button "View Slide Image Biospecimen"
* Select "Details" on the Image Viewer page
* Verify details fields and values
  |field_name                       |value                                  |
  |---------------------------------|---------------------------------------|
  |File_id                          |1fa7d4b5-147a-47a5-8605-3baa0d7be899   |
  |Submitter_id                     |TCGA-GV-A3QI-01A-01-TSA                |
  |Slide_id                         |88217dc9-06a0-4839-8269-85ce0798ef89   |

## Navigate back to Case Summary Page: TCGA-GV-A3QI
* Quick search for "06f936e9-5a90-40d3-b91a-713f2b4e6e11" and go to its page
## Biospecimen Table - Analyte
* Enter "873bc630-7061-4dee-9973-4a4f1b979eaf" in the text box "Biospecimen Search Bar"
* Verify the table "Selection Information Biospecimen" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Analyte ID                             |
    |Analyte UUID                           |
    |Analyte Type                           |
    |Well Number                            |
    |Amount                                 |
    |A260 A280 Ratio                        |
    |Concentration                          |
    |Spectrophotometer Method               |
    |Aliquots                               |
    |TCGA-GV-A3QI-01A-11W                   |
    |873bc630-7061-4dee-9973-4a4f1b979eaf   |
    |Repli-G (Qiagen) DNA                   |

## Biospecimen Table - Aliquot
* Enter "59c24b21-b1e8-4752-bd0a-b48dd802643b" in the text box "Biospecimen Search Bar"
* Verify the table "Selection Information Biospecimen" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Aliquot ID                             |
    |Aliquot UUID                           |
    |Source Center                          |
    |Amount                                 |
    |Concentration                          |
    |TCGA-GV-A3QI-10A-01D-A21Z-08           |
    |59c24b21-b1e8-4752-bd0a-b48dd802643b   |
    |23                                     |
    |0.07                                   |

## Biospecimen Table - TSV Download
* Download "TSV" from "Case Summary Biospecimen Table"
* Read file content from compressed "TSV from Case Summary Biospecimen Table"
* Verify that "TSV from Case Summary Biospecimen Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |project_id                             |
    |sample_id                              |
    |portion_id                             |
    |analytes.submitter_id                  |
    |aliquots.aliquot_id                    |
    |TCGA-GV-A3QI-10A-01W-A226-08           |
    |06f936e9-5a90-40d3-b91a-713f2b4e6e11   |
    |45f1a50a-8d05-48c5-8f82-ba4955117d1b   |
    |TCGA-GV-A3QI-01A-11                    |
    |ffc59e86-1afb-4cdb-8651-191b21ee7876   |
    |TCGA-GV-A3QI-01A-11D-A21Z-08           |
    |analytes.experimental_protocol_type    |
    |analytes.spectrophotometer_method      |
    |mirVana (Allprep DNA) RNA              |
    |analytes.analyte_type_id               |
    |slide_id                               |
    |88217dc9-06a0-4839-8269-85ce0798ef89   |
    |slides.percent_inflam_infiltration     |
    |slides.percent_tumor_nuclei            |
    |portions.submitter_id                  |
    |portions.creation_datetime             |

## Biospecimen Table - Download JSON
* Download "JSON" from "Case Summary Biospecimen Table"
* Read from "JSON from Case Summary Biospecimen Table"
* Verify that "JSON from Case Summary Biospecimen Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |06f936e9-5a90-40d3-b91a-713f2b4e6e11   |
    |TCGA-BLCA                              |
    |TCGA-GV-A3QI                           |
    |074df1cb-b6b7-4d3f-9c0a-fe523b1854e0   |
    |Blood Derived Normal                   |
    |9d5d3fb7-22e5-4ad5-9036-3074e8ac41a8   |
    |3792A65B-C233-4DDF-88F8-8829B73616FE   |
    |873bc630-7061-4dee-9973-4a4f1b979eaf   |
    |32a26490-4b9f-494c-b298-1bd7367aa070   |
    |88217dc9-06a0-4839-8269-85ce0798ef89   |
    |7375d70d-f43f-5589-931e-59e50d743f33   |

## Biospecimen Table - Validate JSON File Fields
  |field_name                               |
  |-----------------------------------------|
  |case_id	                                |
  |project.project_id                       |
  |submitter_id	                            |
  |samples.sample_type_id	                  |
  |samples.portions.portion_id              |
  |samples.portions.analytes.analyte_id     |
  |samples.pathology_report_uuid	          |
  |samples.portions.slides.slide_id         |
  |samples.portions.analytes.aliquots.aliquot_quantity|
* Verify that the "JSON from Case Summary Biospecimen Table" has <field_name> for each object

## Files Table - Count
* Collect "File" Count from Case Summary Header for comparison
* Collect table "Files Case Summary" Item Count for comparison
* Verify "File Count from Case Summary Header" and "Files Case Summary Item Count" are "Equal"

## Files Table - TSV Download
* Search "tcga-blca" in the files table on the case summary page
* Download "TSV" from "Case Summary Files Table"
* Read from "TSV from Case Summary Files Table"
* Verify that "TSV from Case Summary Files Table" has expected information
  |required_info                            |
  |-----------------------------------------|
  |Access                                   |
  |File Name                                |
  |Data Category                            |
  |Data Format                              |
  |Experimental Strategy                    |
  |File Size                                |
  |TCGA-BLCA.6f417a24-f868-4cd1-a60b-37df5f098c76.gene_level_copy_number.v36.tsv|
  |Copy Number Variation                    |
  |BEDPE                                    |
  |RNA-Seq                                  |
  |39925                                    |
  |Controlled                               |
* Verify that "TSV from Case Summary Files Table" does not contain specified information
  |required_info                                                            |
  |-------------------------------------------------------------------------|
  |nationwidechildrens.org_ssf_normal_controls_blca.txt                     |
  |TCGA-GV-A3QI-01A-01-TSA.88217DC9-06A0-4839-8269-85CE0798EF89.svs         |
  |fadf6ebc-28d7-44a5-af0d-518f36cad472.mirbase21.mirnas.quantification.txt |

## Files Table - Download JSON
* Download "JSON" from "Case Summary Files Table"
* Read from "JSON from Case Summary Files Table"
* Verify that "JSON from Case Summary Files Table" has expected information
    |required_info                          |
    |---------------------------------------|
    |00980fb9-5519-4245-bd47-f1858dd69aaf   |
    |Structural Variation                   |
    |06f936e9-5a90-40d3-b91a-713f2b4e6e11   |
    |Genotyping Array                       |
    |Affymetrix SNP 6.0                     |
    |Allele-specific Copy Number Segment    |
    |a1098528-bd94-4d24-bd18-b811fdddd951   |
    |TCGA-BLCA                              |
* Verify that "JSON from Case Summary Files Table" does not contain specified information
  |required_info                                                            |
  |-------------------------------------------------------------------------|
  |nationwidechildrens.org_ssf_normal_controls_blca.txt                     |
  |TCGA-GV-A3QI-01A-01-TSA.88217DC9-06A0-4839-8269-85CE0798EF89.svs         |
  |fadf6ebc-28d7-44a5-af0d-518f36cad472.mirbase21.mirnas.quantification.txt |

## Files Table - Validate JSON File Fields
  |field_name                               |
  |-----------------------------------------|
  |data_format                              |
  |cases.case_id                            |
  |cases.project.project_id                 |
  |access	                                  |
  |file_name                                |
  |file_id                                  |
  |data_type                                |
  |data_category                            |
  |experimental_strategy                    |
  |platform                                 |
  |file_size                                |
* Verify that the "JSON from Case Summary Files Table" has <field_name> for each object

## Files Table - Click Link
* Search "8c239c32-d03a-41fa-9806-2f336cbc6e36" in the files table on the case summary page
* Wait for table "Files Case Summary" body text to appear
    |expected_text                                        |row  |column |
    |-----------------------------------------------------|-----|-------|
    |nationwidechildrens.org_biospecimen.TCGA-GV-A3QI.xml |1    |2      |
* Select value from table "Files Case Summary" by row and column
    |row   |column|
    |------|------|
    |1     |2     |
* Verify the table "File Properties File Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |8c239c32-d03a-41fa-9806-2f336cbc6e36   |
    |nationwidechildrens.org_biospecimen.TCGA-GV-A3QI.xml|
    |3641639dd55ff5eebd20d1ba9286163c       |

## Navigate back again to Case Summary Page: TCGA-GV-A3QI
* Quick search for "06f936e9-5a90-40d3-b91a-713f2b4e6e11" and go to its page

## Most Frequent Somatic Mutations Table - Download TSV
A follow-up will occur in QA-2288 to test filtered download of TSV.
* Download "TSV" from "Mutation Frequency"
* Read from "TSV from Mutation Frequency"
* Verify that "TSV from Mutation Frequency" has expected information
  |required_info                        |
  |-------------------------------------|
  |chr11:g.89047408G>T                  |
  |Non Coding Transcript Exon           |
  |deleterious                          |
  |Substitution                         |
  |THAP12 D449H                         |

## Most Frequent Somatic Mutations Table - Create Set/Cohort, and Click Link
* In table "Most Frequent Somatic Mutations", search the table for "TTN"
* Select value from table "Most Frequent Somatic Mutations" by row and column
  |row   |column|
  |------|------|
  |1     |1     |
  |2     |1     |
  |3     |1     |
  |4     |1     |
  |5     |1     |
  |6     |1     |
* Select "Save/Edit Mutation Set"
* Select "Save as new mutation set" from dropdown menu
* Enter "Case Summary TTN Mutation Set" in the text box "Set Name"
* Select button "Save"
* Is text "Set has been saved." present on the page
* Pause "3" seconds
* In table "Most Frequent Somatic Mutations", search the table for "chr2:g.119882146C>G"
* Wait for table "Most Frequent Somatic Mutations" body text to appear
  |expected_text        |row  |column |
  |---------------------|-----|-------|
  |chr2:g.119882146C>G  |1    |2      |
* Collect button labels in table "Most Frequent Somatic Mutations" for comparison
  |button_label                                     |row  |column |
  |-------------------------------------------------|-----|-------|
  |chr2:g.119882146C>G Affected Cases in TCGA-BLCA  |1    |6      |
* Select value from table "Most Frequent Somatic Mutations" by row and column
  |row   |column|
  |------|------|
  |1     |6     |
* Name the cohort "chr2:g.119882146C>G Affected Cases in TCGA-BLCA" in the Cohort Bar section
* Perform action and validate modal text
  |Action to Perform|Text to validate in modal                                          |Keep or Remove Modal|
  |-----------------|-------------------------------------------------------------------|--------------------|
  |Save             |chr2:g.119882146C>G Affected Cases in TCGA-BLCA has been saved     |Remove Modal        |
* Select value from table "Most Frequent Somatic Mutations" by row and column
  |row   |column|
  |------|------|
  |1     |2     |
* Verify the table "Summary Mutation Summary" is displaying this information
  |text_to_validate                     |
  |------------------------------------ |
  |4ecb8159-bcbf-5900-a1d5-7c94252f6a5a |
  |chr2:g.119882146C>G                  |
  |Single base substitution             |
* Navigate to "Analysis" from "Header" "section"
* Switch cohort to "chr2:g.119882146C>G Affected Cases in TCGA-BLCA" from the Cohort Bar dropdown list
* "chr2:g.119882146C>G Affected Cases in TCGA-BLCA" should be the active cohort
* Collect Cohort Bar Case Count for comparison
* Verify "Cohort Bar Case Count" and "chr2:g.119882146C>G Affected Cases in TCGA-BLCA" are "Equal"
* Navigate to "Manage Sets" from "Header" "section"
* Select item list for set "Case Summary TTN Mutation Set" on Manage Sets page
* Verify the table "Set Information" is displaying this information
  |text_in_table_to_check               |
  |-------------------------------------|
  |10101b20-6578-548d-8b2b-56d311a9fe3b |
  |2e539afe-524a-53fb-92a2-a286357ae3de |
  |0d12a60f-37d8-5721-940d-f78ffff8a8de |
  |5c10528b-51d9-559c-81af-5b0f2097ad82 |
  |0dc64877-8bfa-5000-adc6-fa60cee382c5 |
  |d3a4d7c4-148a-5fc3-95e6-7645271b5281 |
* Close set panel


## Navigate to Case Summary Page: TCGA-EO-A3KX
* Quick search for "fe2e89f7-8f4d-420a-a551-4877cf0fd1d3" and go to its page

## Annotations Table
* In table "Annotations Case Summary", search the table for "064a8d49-9a7a-4667-a757-be5b6d53076d"
* Verify the table "Annotations Case Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |UUID                                   |
    |Entity Type                            |
    |Entity ID                              |
    |Category                               |
    |Classification                         |
    |Created Datetime                       |
    |064a8d49-9a7a-4667-a757-be5b6d53076d   |
    |masked_somatic_mutation                |
    |a4946bbc-5a04-4613-a151-ef9f834b02c0	  |
    |General                                |
    |Redaction                              |
    |2023-03-03T10:25:36.534006-06:00       |

## Annotations Table - Download TSV
* Download "TSV" from "Case Summary Annotations Table"
* Read from "TSV from Case Summary Annotations Table"
* Verify that "TSV from Case Summary Annotations Table" has expected information
    |text_to_validate                       |
    |---------------------------------------|
    |UUID                                   |
    |Entity Type                            |
    |Entity ID                              |
    |Category                               |
    |Classification                         |
    |Created Datetime                       |
    |064a8d49-9a7a-4667-a757-be5b6d53076d   |
    |masked_somatic_mutation                |
    |a4946bbc-5a04-4613-a151-ef9f834b02c0	  |
    |General                                |
    |Redaction                              |
    |2023-03-03T10:25:36.534006-06:00       |
* Verify that "TSV from Case Summary Annotations Table" does not contain specified information
    |text_to_validate                       |
    |---------------------------------------|
    |12dde7d2-fa6c-425b-8f31-1de00f545607   |
    |aggregated_somatic_mutation            |
    |09dc29c3-0541-4563-9861-16ce01f3fd23   |
    |7c54d12d-76f4-4d65-a822-06c67ca6b352   |

## Annotations Table - Download JSON
* Download "JSON" from "Case Summary Annotations Table"
* Read from "JSON from Case Summary Annotations Table"
* Verify that "JSON from Case Summary Annotations Table" has expected information
    |text_to_validate                       |
    |---------------------------------------|
    |064a8d49-9a7a-4667-a757-be5b6d53076d   |
    |a4946bbc-5a04-4613-a151-ef9f834b02c0   |
    |DEV-1606, redaction on parent node did not get applied to this node|
    |Redaction                              |
    |General                                |
    |2023-03-03T10:25:36.534006-06:00       |
    |TCGA-UCEC                              |
    |TCGA-EO-A3KX                           |
    |0d4b1931-94ff-4dd3-8646-b97eaed63d6d   |
* Verify that "JSON from Case Summary Annotations Table" does not contain specified information
    |text_to_validate                       |
    |---------------------------------------|
    |12dde7d2-fa6c-425b-8f31-1de00f545607   |
    |aggregated_somatic_mutation            |
    |09dc29c3-0541-4563-9861-16ce01f3fd23   |
    |7c54d12d-76f4-4d65-a822-06c67ca6b352   |

## Annotations Table - Validate JSON File Fields
  |field_name                               |
  |-----------------------------------------|
  |annotation_id	                          |
  |entity_submitter_id	                    |
  |notes	                                  |
  |entity_type	                            |
  |case_id	                                |
  |project.project_id                       |
  |project.program.name                     |
  |classification	                          |
  |entity_id	                              |
  |category	                                |
  |created_datetime                         |
  |status                                   |
  |case_submitter_id	                      |
* Verify that the "JSON from Case Summary Annotations Table" has <field_name> for each object

## Annotations Table - Click Link
* Select value from table "Annotations Case Summary" by row and column
    |row   |column|
    |------|------|
    |1     |1     |
* Verify the table "Summary Annotation Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |Annotation UUID                        |
    |064a8d49-9a7a-4667-a757-be5b6d53076d   |
* Verify the table "Notes Annotation Summary" is displaying this information
    |text_to_validate                       |
    |---------------------------------------|
    |DEV-1606, redaction on parent node did not get applied to this node|
