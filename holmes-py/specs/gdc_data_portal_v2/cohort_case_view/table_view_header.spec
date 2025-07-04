# Cohort Case View - Table View Header
Date Created        : 10/07/2024
Version			    : 1.0
Owner		        : GDC QA
Description		    : Validate Cohort Header
Test-Case           : PEAR-1052

tags: gdc-data-portal-v2, regression, cohort-bar, case-view, clinical-biospecimen-download

## Navigate to Table View
* On GDC Data Portal V2 app
* Navigate to "Analysis" from "Header" "section"
* Expand or collapse the cohort bar
* Go to tab "Table View" in Cohort Case View
* Unpin the Cohort Bar

## Biospecimen - TSV
* In table "Cases", search the table for "MMRF_11"
* Download "TSV" from "Cohort Table View Biospecimen"
* Read file content from compressed "TSV from Cohort Table View Biospecimen"
* Verify that "TSV from Cohort Table View Biospecimen" has expected information
    |required_info                          |
    |---------------------------------------|
    |aliquots.aliquot_id                    |
    |aliquots.aliquot_volume                |
    |no_matched_normal_wxs                  |
    |aliquots.selected_normal_targeted_sequencing|
    |7c0feb4c-3ce1-529a-98a6-3e48e1d06a8a   |
    |ba678f26-a87d-4da0-bcee-d3a6caffae53   |
    |MMRF_1152_1_BM_CD138pos_T2_TSE61_K03194|
    |analyte_quantity                       |
    |ribosomal_rna_28s_16s_ratio            |
    |well_number                            |
    |172f1884-c9d3-529b-a885-b5a37442a7e3   |
    |de08e2be-078b-58b3-90a6-7546145af589   |
    |slide_id                               |
    |slides.submitter_id                    |
    |percent_neutrophil_infiltration        |
    |prostatic_chips_total_count            |
    |portion_id                             |
    |portions.submitter_id                  |
    |portions.portion_number                |
    |8742f2f2-fe95-5f62-9ef1-8241488cd623   |
    |331d5f27-bfe8-58f6-98b9-f54376f15c96   |
    |sample_id                              |
    |samples.submitter_id                   |
    |days_to_sample_procurement             |
    |intermediate_dimension                 |
    |samples.preservation_method            |
    |samples.time_between_clamping_and_freezing|
    |samples.tumor_code_id                  |
    |ccff20ac-6979-4602-ae3e-bb8fb3968a61   |
    |-27                                    |
    |Primary Blood Derived Cancer - Bone Marrow|
    |Bone Marrow NOS                        |
    |MMRF_1195_1_PB_Whole                   |
* Verify that "TSV from Cohort Table View Biospecimen" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |TCGA                                   |

## Biospecimen - JSON
* In table "Cases", search the table for "MMRF_1186"
* Download "JSON" from "Cohort Table View Biospecimen"
* Read from "JSON from Cohort Table View Biospecimen"
* Verify that "JSON from Cohort Table View Biospecimen" has expected information
    |required_info                          |
    |---------------------------------------|
    |e04417a0-1644-410d-b771-d5d097500a49   |
    |MMRF_1186                              |
    |MMRF-COMMPASS                          |
    |0ee3be21-7250-4f1c-a5c0-6b2ba89b2f4a   |
    |MMRF_1186_2_BM_CD138pos                |
    |c15574aa-1229-5e79-97c8-e3f0052c1442   |
    |80b55c55-e704-439a-919f-9661cc94ba8d   |
    |MMRF_1186_2_BM_CD138pos_T1_KHS5U_L14411|
    |ff9daa57-3c56-5a00-b2d7-6fb061caf3fa   |
    |2018-01-01T12:08:49.642883-06:00       |
    |837e028e-39f2-51f2-9ae7-6242bb0dff0a   |
* Verify that "JSON from Cohort Table View Biospecimen" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |TCGA                                   |
    |MMRF_1108                              |

## Biospecimen - Validate JSON Fields
  |field_name                                       |
  |-------------------------------------------------|
  |case_id                                          |
  |project.project_id                               |
  |submitter_id                                     |
  |samples.tumor_descriptor	                        |
  |samples.preservation_method                      |
  |samples.sample_id                                |
  |samples.portions.analytes.analyte_id	            |
  |samples.portions.analytes.aliquots.aliquot_id	|
  |samples.portions.portion_id	                    |
* Verify that the "JSON from Cohort Table View Biospecimen" has <field_name> for each object

## Clinical - TSV
* In table "Cases", search the table for "MMRF_28"
* Download "TSV" from "Cohort Table View Clinical"
* Read file content from compressed "TSV from Cohort Table View Clinical"
* Verify that "TSV from Cohort Table View Clinical" has expected information
    |required_info                          |
    |---------------------------------------|
    |follow_ups.days_to_follow_up           |
    |follow_ups.evidence_of_recurrence_type |
    |follow_ups.submitter_id                |
    |follow_ups.pregnancy_outcome           |
    |molecular_tests.laboratory_test        |
    |f7253ebc-424c-4506-a5f4-727c819d50fd   |
    |d185b62d-3811-43c1-8b50-0823e06646fa   |
    |MMRF_2838_molecular_test6              |
    |Test Value Reported                    |
    |30.79999                               |
    |pathology_details.lymphatic_invasion_present|
    |pathology_details.submitter_id         |
    |pathology_details.size_extraocular_nodule|
    |pathology_details.zone_of_origin_prostate|
    |exposures.alcohol_drinks_per_day       |
    |demographic.occupation_duration_years  |
    |family_histories.relationship_primary_diagnosis|
    |family_histories.relatives_with_cancer_history_count|
    |demographic.weeks_gestation_at_birth   |
    |diagnoses.child_pugh_classification    |
    |diagnoses.first_symptom_longest_duration|
    |treatments.prescribed_dose_units       |
    |1d485291-469a-4d6e-950e-06b974f4445d   |
    |MMRF_2832                              |
    |9732/3                                 |
    |First line of therapy                  |
    |Stem Cell Transplantation, Autologous  |
    |-16767                                 |
    |MMRF-COMMPASS                          |
    |cases.case_id                          |
    |cases.consent_type                     |
    |cases.days_to_consent                  |
    |cases.days_to_lost_to_followup         |
    |cases.disease_type                     |
    |cases.index_date                       |
    |cases.lost_to_followup                 |
    |cases.primary_site                     |
    |cases.submitter_id                     |

* Verify that "TSV from Cohort Table View Clinical" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |TCGA                                   |
    |MMRF_1108                              |
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

## Clinical - JSON
* In table "Cases", search the table for "MMRF_187"
* Download "JSON" from "Cohort Table View Clinical"
* Read from "JSON from Cohort Table View Clinical"
* Verify that "JSON from Cohort Table View Clinical" has expected information
    |required_info                          |
    |---------------------------------------|
    |Plasma Cell Tumors                     |
    |Hematopoietic and reticuloendothelial systems|
    |First Treatment                        |
    |cd63c5c4-308e-4d7b-bd41-22946218fcef   |
    |Melanoma                               |
    |8e39bdf0-044d-4e25-a7ac-14239117e580   |
    |0a8eb97e-b621-4bae-97c1-832462fa7e32   |
    |MMRF_1876_molecular_test78             |
    |Platelets                              |
    |MMRF_1876_followup5                    |
    |2018-07-19T16:28:52.616328-05:00       |
    |MMRF-COMMPASS                          |
    |MMRF_1879                              |
    |MMRF_1879_diagnosis1                   |
    |5e7224a2-8639-45e8-b470-5d034d2483dc   |
    |First line of therapy                  |
    |Bortezomib                             |
    |0b3b060b-122c-4a75-8c1d-2ac374e65aec   |
    |black or african american              |
    |c04a3091-378f-47a1-8150-a804519d0031   |
    |323                                    |
* Verify that "JSON from Cohort Table View Clinical" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |TCGA                                   |
    |MMRF_1108                              |

## Clinical - TSV - Verify Additional Case Properties
* In table "Cases", search the table for "01b7ac35-49e3-48db-a081-2971b807445f"
* Download "TSV" from "Cohort Table View Clinical"
* Read file content from compressed "TSV from Cohort Table View Clinical"
* Verify that "TSV from Cohort Table View Clinical" has expected information
    |required_info                          |
    |---------------------------------------|
    |Informed Consent                       |
    |-3                                     |
    |702                                    |
    |Ductal and Lobular Neoplasms           |
    |Diagnosis                              |
    |Yes                                    |
    |Pancreas                               |
    |C3N-02696                              |

## Clinical - JSON - Verify Additional Case Properties
* Download "JSON" from "Cohort Table View Clinical"
* Read from "JSON from Cohort Table View Clinical"
* Verify that "JSON from Cohort Table View Clinical" has expected information
    |required_info                          |
    |---------------------------------------|
    |Informed Consent                       |
    |-3                                     |
    |702                                    |
    |Ductal and Lobular Neoplasms           |
    |Diagnosis                              |
    |Yes                                    |
    |Pancreas                               |
    |C3N-02696                              |

## Clinical - Validate JSON Fields
  |field_name                                       |
  |-------------------------------------------------|
  |case_id                                          |
  |consent_type                                     |
  |days_to_consent                                  |
  |days_to_lost_to_followup                         |
  |disease_type                                     |
  |index_date                                       |
  |lost_to_followup                                 |
  |primary_site                                     |
  |submitter_id                                     |
* Verify that the "JSON from Cohort Table View Clinical" has <field_name> for each object

## Cohort Table - TSV
* In table "Cases", search the table for "TCGA-EK"
* Download "TSV" from "Cohort Table View"
* Read from "TSV from Cohort Table View"
* Verify that "TSV from Cohort Table View" has expected information
    |required_info                          |
    |---------------------------------------|
    |demographic.days_to_death              |
    |diagnoses.0.age_at_diagnosis           |
    |disease_type                           |
    |submitter_id                           |
    |primary_site                           |
    |summary.experimental_strategies.0.experimental_strategy|
    |summary.file_count                     |
    |16363ef7-899d-4927-a74b-e7c03b3d8af1   |
    |d4510801-8ce4-4f10-8b31-25c0a39b5135   |
    |Squamous cell carcinoma, large cell, nonkeratinizing, NOS|
    |18630                                  |
    |f92a994d-b080-42c9-966f-34c1000865e3   |
    |TCGA-CESC                              |
    |TCGA                                   |
    |TCGA-EK-A3GN                           |
    |Genotyping Array                       |
    |Reverse Phase Protein Array            |
    |Tissue Slide                           |
    |63                                     |
* Verify that "TSV from Cohort Table View" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |MMRF                                   |


## Cohort Table - JSON
* In table "Cases", search the table for "MMRF_11"
* Download "JSON" from "Cohort Table View"
* Read from "JSON from Cohort Table View"
* Verify that "JSON from Cohort Table View" has expected information
    |required_info                          |
    |---------------------------------------|
    |RNA-Seq                                |
    |WGS                                    |
    |WXS                                    |
    |Hematopoietic and reticuloendothelial systems|
    |Plasma Cell Tumors                     |
    |8bb96689-fe0a-45b5-af07-3a6c30480927   |
    |MMRF_1171                              |
    |MMRF-COMMPASS                          |
    |MMRF                                   |
    |29995                                  |
    |white                                  |
    |male                                   |
    |not hispanic or latino                 |
    |Dead                                   |
    |female                                 |
    |372                                    |
    |Alive                                  |
    |6356f1f3-7d56-4464-b70c-84eb69b7262c   |
* Verify that "JSON from Cohort Table View" does not contain specified information
    |required_info                          |
    |---------------------------------------|
    |FM-AD                                  |
    |TARGET-AML                             |
    |MATCH                                  |
    |TCGA                                   |

## Cohort Table - Validate JSON Fields
  |field_name                                       |
  |-------------------------------------------------|
  |summary.file_count                               |
  |primary_site	                                    |
  |disease_type	                                    |
  |case_id	                                        |
  |submitter_id	                                    |
  |project.project_id                               |
  |project.program.name                             |
  |diagnoses.age_at_diagnosis                       |
  |demographic.race                                 |
  |demographic.gender                               |
  |demographic.ethnicity                            |
  |demographic.vital_status                         |
  |demographic.days_to_death                        |
* Verify that the "JSON from Cohort Table View" has <field_name> for each object
