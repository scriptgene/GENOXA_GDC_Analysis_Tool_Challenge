# Cohort Bar - Save, Switch, Delete Cohort
Date Created    : 03/29/2023
Version			    : 1.0
Owner		        : GDC QA
Description		  : Test Cohort Bar - saving, switching, and deleting a cohort
Test-case       : PEAR-492,PEAR-493,PEAR-498

tags: gdc-data-portal-v2, regression, cohort-bar,smoke-test

## Navigate to Cohort Builder
* On GDC Data Portal V2 app
* Navigate to "Cohort" from "Header" "section"

## Save a cohort for the first time
* Select "Add" from the Cohort Bar
* Is modal with text "Unsaved_Cohort has been created" present on the page and "Remove Modal"
* Make the following selections from "Demographic" tab on the Cohort Builder page
  |facet_name       |selection            |
  |-----------------|---------------------|
  |Gender           |female               |
  |Vital Status     |alive                |
* Is text "Cohort not saved" present on the Cohort Bar
* Select "Save" from the Cohort Bar
* Name the cohort "test save 1" in the Cohort Bar section
* Perform action and validate modal text
  |Action to Perform|Text to validate in modal                  |Keep or Remove Modal|
  |-----------------|-------------------------------------------|--------------------|
  |Save             |Cohort has been saved                      |Remove Modal        |
* "test save 1" should be the active cohort
* Is text "Cohort not saved" not present on the Cohort Bar

## Save a cohort for the second time
* Make the following selections from "General" tab on the Cohort Builder page
  |facet_name       |selection                    |
  |-----------------|-----------------------------|
  |Disease Type     |adenomas and adenocarcinomas |
* Is text "Changes not saved" present on the Cohort Bar
* Select "Save" from the Cohort Bar
* Is text "You cannot undo this action." present on the page
* Perform action and validate modal text
  |Action to Perform|Text to validate in modal                  |Keep or Remove Modal|
  |-----------------|-------------------------------------------|--------------------|
  |Save             |Cohort has been saved                      |Remove Modal        |
* "test save 1" should be the active cohort
* Is text "Changes not saved" not present on the Cohort Bar

## Switch cohorts
* Select "Add" from the Cohort Bar
* Is modal with text "Unsaved_Cohort has been created" present on the page and "Remove Modal"
* Switch cohort to "test save 1" from the Cohort Bar dropdown list
* "test save 1" should be the active cohort
* Is text "Cohort not saved" not present on the Cohort Bar
* Switch cohort to "Unsaved_Cohort" from the Cohort Bar dropdown list
* "Unsaved_Cohort" should be the active cohort
* Is text "Cohort not saved" present on the Cohort Bar
* Switch cohort to "test save 1" from the Cohort Bar dropdown list
* "test save 1" should be the active cohort

## Delete a saved cohort
* Select "Delete" from the Cohort Bar
* Perform action and validate modal text
  |Action to Perform|Text to validate in modal                  |Keep or Remove Modal|
  |-----------------|-------------------------------------------|--------------------|
  |Delete           |test save 1                                |Remove Modal        |
* Switch cohort to "Unsaved_Cohort" from the Cohort Bar dropdown list
* "Unsaved_Cohort" should be the active cohort

## Delete an unsaved cohort
* Expand or contract a facet from "Exposure" tab on the Cohort Builder page
  |facet_name               |selection                    |
  |-------------------------|-----------------------------|
  |Tobacco Smoking Status   |plus-icon                    |
* Make the following selections from "Exposure" tab on the Cohort Builder page
  |facet_name               |selection                    |
  |-------------------------|-----------------------------|
  |Tobacco Smoking Status   |lifelong non-smoker          |
* Is text "Cohort not saved" present on the Cohort Bar
* Select "Delete" from the Cohort Bar
* Perform action and validate modal text
  |Action to Perform|Text to validate in modal                  |Keep or Remove Modal|
  |-----------------|-------------------------------------------|--------------------|
  |Delete           |has been deleted                           |Remove Modal        |

## Random filters cohort creation
* Create and save cohorts with randomly assigned filters
  |cohort_name              |number_of_filters|
  |-------------------------|-----------------|
  |random_filters_cohort_1  |2                |

## Shortened cohort saving
* Create and save a cohort named "quick_cohort_step_1" with these filters
  |tab_name               |facet_name           |selection            |
  |-----------------------|---------------------|---------------------|
  |Demographic            |Gender               |female               |
  |General Diagnosis      |Ajcc Clinical Stage  |stage ivc            |
  |General                |Program              |TCGA                 |
