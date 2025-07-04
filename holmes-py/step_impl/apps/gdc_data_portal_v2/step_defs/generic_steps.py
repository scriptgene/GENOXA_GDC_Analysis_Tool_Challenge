import json
import re
import tarfile
import time

from datetime import datetime as dt
from uuid import UUID, uuid4

from getgauge.python import step, before_spec, after_spec, before_suite, data_store

from ..app import GDCDataPortalV2App
from ....base.webdriver import WebDriver
from ....base.utility import Utility


@step("Pause <sleep_time> seconds")
def pause_10_seconds(sleep_time):
    time.sleep(int(sleep_time))


@before_suite
def start_app():
    global APP
    APP = GDCDataPortalV2App(WebDriver.page)


@before_suite
def navigate_to_app():
    APP.navigate()
    APP.modal.accept_warning()


@before_suite
def setup_test_run():
    """
    Before the start of tests this function will run once. When accessing the data portal for the
    first time, an unsaved cohort will be created. There can only be one unsaved cohort in the
    data portal at a time. That means the user cannot create another cohort (in the cohort bar), until
    the initial one has been saved. This function saves that initial cohort so subsequent tests can run, as some
    depend on a starting condition that they be able to create a cohort.
    """
    APP.analysis_center_page.visit()
    APP.header_section.wait_for_page_to_load("analysis")
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    time.sleep(2)
    APP.cohort_bar.click_cohort_bar_button("Save")
    APP.shared.click_text_option_from_dropdown_menu("Save")
    APP.shared.send_text_into_text_box("never_use_this_cohort_name", "Name Input Field")
    APP.shared.click_button_in_modal_with_displayed_text_name("Save")
    APP.cohort_bar.wait_for_text_in_temporary_message(
        "Cohort has been saved", "Remove Modal"
    )
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    time.sleep(2)


@before_spec("<controlled-access>")
def before_spec_hook():
    """
    Before each spec is ran, automation will check if the file contains the tag "controlled-access".
    If it does, automation will check the data portal to see if the user is logged in. If not, automation
    will click the login button and we will manually login to the data portal.

    Gauge does not allow for this kind of tag checking in before_suite, so this is my workaround. To do
    it before spec in an 'if' statement.
    """
    APP.header_section.login_to_data_portal_if_possible()


@after_spec
def setup_next_spec_run():
    """
    After each spec file's execution this function will run. The intention is to
    clear the active cohort filters and setup the next spec run.

    First, we go to the analysis center. If a test found a bug in the data portal the next test
    may not execute correctly, because the previous test ended in an unexpected and difficult to get
    out of place. Then, we check to see if there is the 'No filters currently applied' text is present.
    If not, we click the 'Clear All' button to remove the active cohort filters.
    Finally, we wait to see the the text confirming there are no active cohort filters present.
    """
    APP.analysis_center_page.visit()
    APP.header_section.wait_for_page_to_load("analysis")
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    if not APP.shared.is_no_active_cohort_filter_text_present():
        APP.shared.clear_active_cohort_filters()


@step("On GDC Data Portal V2 app")
def navigate_to_app():
    APP.navigate()
    APP.modal.accept_warning()


@step("Go to <page_name> page")
def go_to_page(page_name):
    pages = {"Analysis": APP.analysis_center_page.visit()}
    pages[page_name]


@step("Navigate to <target> from <source> <target_type>")  # will have click actions
def navigate_to_page_in_page(target, source, target_type):
    sources = {
        "Header": {"section": APP.header_section.navigate_to_main_pages},
        "Home Page": {"app": APP.home_page.navigate_to_app},
        "Repository": {"app": APP.repository_page.click_button},
        "Analysis": {"app": APP.analysis_center_page.navigate_to_app},
        "Cohort Builder": {"app": APP.cohort_builder_page.click_button},
    }
    sources.get(source)[target_type](target)


@step("Verify that the <text> text is displayed on <source> <target_type>")
def verify_text_on_page(text, source, target_type):
    sources = {
        "Repository": {"app": APP.repository_page.get_title},
        "Add a Custom Filter": {
            "modal": APP.repository_page.get_text_on_add_custom_filter_modal
        },
    }
    text_value = sources.get(source)[target_type](text)
    assert (
        text == text_value
    ), f"Unexpected title detected: looking for {text}, but got {text_value}"


@step("Verify <statistic_1> and <statistic_2> are <equal_or_not_equal>")
def verify_compared_statistics_are_equal_or_not_equal(
    statistic_1, statistic_2, equal_or_not_equal
):
    """
    verify_compared_statistics_are_equal_or_not_equal compares two previously stored statistics to one another.
    Asserts if they are equal or not equal based on spec file input.

    :param statistic_1: The first statistic to compare. It is the name it is stored under in data_store.spec
    :param statistic_2: The second statistic to compare. It is the name it is stored under in data_store.spec
    :param equal_or_not_equal: If the compared statistics should be equal or not
    :return: N/A
    """
    # Get first statistic to compare
    first_statistic_string = data_store.spec[f"{statistic_1}"]
    first_statistic_string = APP.shared.strip_string_for_comparison(
        first_statistic_string
    )

    # Get second statistic to compare
    second_statistic_string = data_store.spec[f"{statistic_2}"]
    second_statistic_string = APP.shared.strip_string_for_comparison(
        second_statistic_string
    )

    equal_or_not_equal = equal_or_not_equal.lower()
    if equal_or_not_equal == "equal":
        assert (
            first_statistic_string == second_statistic_string
        ), f"The first statistic {statistic_1}'s value '{first_statistic_string}' and second statistic {statistic_2}'s value '{second_statistic_string}' does NOT match"
    elif equal_or_not_equal == "not equal":
        assert (
            first_statistic_string != second_statistic_string
        ), f"The first statistic {statistic_1}'s value '{first_statistic_string}' and second statistic {statistic_2}'s value '{second_statistic_string}' does match when it should NOT"


@step("Close the modal")
def close_the_modal():
    APP.shared.click_close_modal_button()
    time.sleep(0.5)

@step("Close the message")
def close_the_modal():
    """Clicks the 'x' to close temporary message"""
    APP.shared.click_close_temporary_message()

@step("Download <file> from <source>")
def download_file_at_file_table(file: str, source: str):
    sources = {
        "Browse Annotations": APP.browse_annotations.click_annotation_download_button,
        "Cart Items": APP.shared.click_button_data_testid_normalize,
        "Cart Header": APP.shared.click_button_with_displayed_text_name,
        "Cart Header Dropdown": APP.shared.click_text_option_from_dropdown_menu,
        "Case Summary Annotations Table": APP.case_summary_page.click_annotations_table_download_json_or_tsv_button,
        "Case Summary Biospecimen Table": APP.case_summary_page.click_biospecimen_table_download_button,
        "Case Summary Clinical Table": APP.case_summary_page.click_clinical_table_download_button,
        "Case Summary Files Table": APP.case_summary_page.click_files_table_download_json_or_tsv_button,
        "Case Summary Files Table File": APP.case_summary_page.click_files_table_download_file_button,
        "CDAVE Analysis Card TSV": APP.clinical_data_analysis.click_tsv_button_on_analysis_card,
        "Cohort Bar": APP.cohort_bar.click_cohort_bar_button,
        "Cohort Case View Files": APP.cohort_case_view_page.click_files_and_dropdown_option_cases_view,
        "Cohort Comparison": APP.cohort_comparison_page.click_download_tsv_button_on_analysis_card_cohort_comparison,
        "Cohort Table View": APP.cohort_case_view_page.click_table_view_button,
        "Cohort Summary View Biospecimen": APP.cohort_case_view_page.click_biospecimen_summary_view,
        "Cohort Summary View Clinical": APP.cohort_case_view_page.click_clinical_summary_view,
        "Cohort Table View Biospecimen": APP.cohort_case_view_page.click_biospecimen_table_view,
        "Cohort Table View Clinical": APP.cohort_case_view_page.click_clinical_table_view,
        "Controlled Access Modal": APP.modal.click_download_button,
        "File Summary": APP.file_summary_page.click_download_button,
        "File Summary Annotation Table": APP.file_summary_page.click_annotation_table_download_option,
        "File Summary File Versions": APP.file_summary_page.click_file_version_download_option,
        "Gene Summary Cancer Distribution": APP.gene_summary_page.click_cancer_distribution_button,
        "Graph Dropdown": APP.shared.click_text_option_from_dropdown_data_image,
        "Manage Sets": APP.manage_sets_page.click_on_download_for_set,
        "Most Frequent Somatic Mutations": APP.shared.click_button_in_most_frequent_somatic_mutations,
        "Mutation Frequency": APP.mutation_frequency_page.click_table_download_button,
        "Mutation Summary Consequences": APP.mutation_summary_page.click_consequences_button,
        "Mutation Summary Cancer Distribution": APP.mutation_summary_page.click_cancer_distribution_button,
        "Page": APP.shared.click_button_with_displayed_text_name,
        "Projects": APP.projects_page.click_button,
        "Project Summary": APP.project_summary_page.click_button,
        "Project Summary Biospecimen": APP.project_summary_page.click_biospecimen_download_button,
        "Project Summary Clinical": APP.project_summary_page.click_clinical_download_button,
        "Project Summary Annotations": APP.project_summary_page.click_annotation_download_button,
        "Repository": APP.repository_page.click_button,
        "Set Operations": APP.set_operations_page.click_download_tsv_button_set_operations,
        "Set Operations Union Row": APP.set_operations_page.click_union_row_download_tsv_button_set_operations,
    }
    driver = WebDriver.page
    with driver.expect_download(timeout=90000) as download_info:
        # Allows using sources without passing in contents of <file> as a parameter
        if file.lower() == "file":
            sources.get(source)()
        else:
            sources.get(source)(file)
    download = download_info.value
    file_path = f"{Utility.parent_dir()}/downloads/{dt.timestamp(dt.now())}_{download.suggested_filename}"
    download.save_as(file_path)
    data_store.spec[f"{file} from {source}"] = file_path

@step("In table <table_name> on row <row_identifier> select <button_id> to download file <file_name_to_store>")
def click_button_in_table_download_file(table_name:str, row_identifier:str, button_id: str, file_name_to_store:str):
    """
    click_button_in_table_download_file In specified table with given identifier click the download button.
    A file downloads and we store the information using given identifier.

    :param table_name: The table ID to access.
    :param row_identifier: A unique piece of text to identify which row we want to click the button.
    :param button_id: The button ID to click.
    :param file_name_to_store: Identifier used to store contents of file.
    """
    driver = WebDriver.page
    with driver.expect_download(timeout=90000) as download_info:
        APP.file_summary_page.click_button_in_table(table_name, row_identifier, button_id)
    download = download_info.value
    file_path = f"{Utility.parent_dir()}/downloads/{dt.timestamp(dt.now())}_{download.suggested_filename}"
    download.save_as(file_path)
    data_store.spec[f"{file_name_to_store}"] = file_path


@step("Upload <file_name> <extension> from <folder_name> in <source> through <button>")
def upload_file(
    file_name: str, extension: str, folder_name: str, source: str, button: str
):
    """
    upload_file performs an action to launch a file explorer,
    and then uploads a specified file located in the resources folder.

    :param file_name: The name of the file being uploaded
    :param extension: The extension of the file being uploaded
    :param folder_name: Name of the folder that contains the file being uploaded.
    :param source: Specifies what function is causing the action to launch the file explorer
    :param button: ID or Name of the button that is being clicked to launch the file explorer
    :return: N/A
    """
    sources = {
        "Cohort Bar Import": APP.cohort_bar.click_import_cohort_browse,
        "Cohort Filter Import": APP.cohort_bar.click_import_cohort_browse,
        "Mutation Frequency Custom Filter": APP.mutation_frequency_page.click_custom_filter_import_browse,
        "Manage Sets Import": APP.manage_sets_page.click_browse_import_set,
    }
    driver = WebDriver.page
    with driver.expect_file_chooser(timeout=60000) as file_chooser_info:
        # Allows using sources without passing in contents of <file> as a parameter
        sources.get(source)(button)
    file_chooser = file_chooser_info.value
    file_name = file_name.lower().replace(" ", "_")
    folder_name = folder_name.lower().replace(" ", "_")
    file_path = (
        f"{Utility.parent_dir()}/resources/{folder_name}/{file_name}.{extension}"
    )
    file_chooser.set_files(file_path)


@step("Read from <file_type>")
def read_from_file(file_type):
    with open(data_store.spec[file_type], "r+") as f:
        data_store.spec[f"{file_type} contents"] = f.read()


@step("Read file content from compressed <file_type>")
def read_file_content_from_compressed_file(file_type):
    """Used for tar.gz files. Typically seen with file or multiple file downloads"""
    tar = tarfile.open(data_store.spec[file_type], "r:gz")
    all_files_content = ""
    for member in tar.getmembers():
        f = tar.extractfile(member)
        single_file_content = f.read()
        all_files_content += str(single_file_content)
    data_store.spec[f"{file_type} contents"] = all_files_content


@step("Read metadata from compressed <file_type>")
def read_metadata_from_compressed_file(file_type):
    """Used for tar.gz files. Typically seen with file or multiple file downloads"""
    tar = tarfile.open(data_store.spec[file_type], "r:gz")
    tar_list = tar.getmembers()
    # The first 'member' is always the metadata file
    tar_list_metadata = tar.extractfile(tar_list[0])
    metadata_content = tar_list_metadata.read()
    data_store.spec[f"{file_type} contents"] = str(metadata_content)


@step("Verify that <file_type> has expected information <table>")
def verify_file_content(file_type, table):
    """Checks if specified information is inside collected content from read-in files"""
    for k, v in enumerate(table):
        assert (
            v[0] in data_store.spec[f"{file_type} contents"]
        ), f"'{v[0]}' is NOT found in the file"


@step("Verify that <file_type> has expected information from collected data <table>")
def verify_file_content(file_type, table):
    """Checks if collected information is inside content from read-in files"""
    for k, v in enumerate(table):
        # Get first statistic to compare
        collected_data_string = data_store.spec[f"{v[0]}"]
        collected_data_string = APP.shared.strip_string_for_comparison(
            collected_data_string
        )
        if collected_data_string == "--":
            assert (
                "0" in data_store.spec[f"{file_type} contents"]
            ), f"Collected info '{v[0]}' with value '0' is NOT found in the file"
        else:
            assert (
                collected_data_string in data_store.spec[f"{file_type} contents"]
            ), f"Collected info '{v[0]}' with value '{collected_data_string}' is NOT found in the file"

@step("Verify that <file_type> has exact expected information from collected data <table>")
def verify_file_content(file_type, table):
    """Checks if collected information is inside content from read-in files. Does not strip string for comparison."""
    for k, v in enumerate(table):
        # Get first statistic to compare
        collected_data_string = data_store.spec[f"{v[0]}"]
        if collected_data_string == "--":
            assert (
                "0" in data_store.spec[f"{file_type} contents"]
            ), f"Collected info '{v[0]}' with value '0' is NOT found in the file"
        else:
            assert (
                collected_data_string in data_store.spec[f"{file_type} contents"]
            ), f"Collected info '{v[0]}' with value '{collected_data_string}' is NOT found in the file"

@step("Verify that <file_type> does not contain specified information <table>")
def verify_content_is_not_in_file(file_type, table):
    """Checks if specified information is NOT inside collected content from read-in files"""
    for k, v in enumerate(table):
        assert (
            v[0] not in data_store.spec[f"{file_type} contents"]
        ), f"'{v[0]}' is found in the file when it's unexpected"


@step("Verify that the <file_type> has <field_name> for each object")
def verify_file_has_expected_field_names(file_type, field_name):
    fails = {}
    if "JSON" in file_type:
        if "annotations.annotation_id" in field_name:
            return
        json_arr = json.loads(data_store.spec[f"{file_type} contents"])
        field_names = field_name.split(".")
        for json_obj in json_arr:
            if len(field_names) > 1:
                Utility.validate_json_key_exists(
                    Utility.flatten_json(json_obj), field_name, fails
                )
            elif "platform" in field_name:
                if (
                    not json_obj["data_format"] == "VCF"
                    and not json_obj["data_format"] == "MAF"
                    and not json_obj["data_format"] == "TSV"
                ):
                    Utility.validate_json_key_exists(json_obj, field_name, fails)
            elif "experimental_strategy" in field_name:
                if (
                    not json_obj["data_category"] == "Clinical"
                    and not json_obj["data_category"] == "Biospecimen"
                ):
                    Utility.validate_json_key_exists(json_obj, field_name, fails)
            else:
                Utility.validate_json_key_exists(json_obj, field_name, fails)
    elif "TSV" in file_type:
        pass
    assert not fails, f"{file_type} validation failed!\nFails: {fails}"


@step("Verify these items are not on the page <table>")
def verify_item_does_not_appear_on_page(table):
    """Verifies if the given data-testid is not present on the page"""
    for k, v in enumerate(table):
        data_testid_to_check = APP.shared.normalize_button_identifier(v[0])

        is_filter_visible = APP.shared.is_data_testid_present(data_testid_to_check)
        assert is_filter_visible == False, f"The item '{v[0]}' is visible when it should NOT be"

        # Also check exactly the way it was given as test data
        is_filter_visible = APP.shared.is_data_testid_present(v[0])
        assert is_filter_visible == False, f"The item '{v[0]}' is visible when it should NOT be"

@step("Verify these custom filters <are_or_are_not> present in the properties table <table>")
def verify_custom_filters_presence(are_or_are_not:str, table):
    """Verifies if the given custom filter is present or not"""
    for k, v in enumerate(table):
        is_custom_filter_visible = APP.shared.is_custom_filter_in_properties_table_present(v[0])
        if are_or_are_not.lower() == "are":
            assert is_custom_filter_visible, f"In properties table, the custom filter '{v[0]}' is NOT visible"
        elif are_or_are_not.lower() == "are not":
            assert is_custom_filter_visible == False, f"In properties table, the custom filter '{v[0]}' is present when it should NOT be"

@step("Verify presence of filter card <table>")
def make_cohort_builder_selections(table):
    for k, v in enumerate(table):
        is_filter_visible = APP.shared.is_filter_card_present(v[0])
        assert is_filter_visible, f"The filter card '{v[0]}' is NOT visible"

@step("Verify expected custom filters <are_or_are_not> present in filter card <table>")
def validate_custom_filter_text_on_facet_card(are_or_are_not:str, table):
    """
    validate_custom_filter_text_on_facet_card validates given text should be present
    or not present on filter card. Used with cards that allow entry of text as a filter
    e.g. Case ID, Mutated Gene, Cases Submitter Id, etc.

    :param are_or_are_not: Check if text should be present or should not be present. Acceptable inputs: "are", "are not"
    :param v[0]: The facet card name to check
    :param v[1]: Custom filter text to check on the facet card
    """
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        is_custom_filter_visible = APP.shared.is_filter_card_custom_filter_text_present(v[0], v[1])
        if are_or_are_not.lower() == "are":
            assert is_custom_filter_visible, f"In filter card '{v[0]}', the custom filter text '{v[1]}' is NOT visible"
        elif are_or_are_not.lower() == "are not":
            assert is_custom_filter_visible == False, f"In filter card '{v[0]}', the custom filter text '{v[1]}' is present when it should NOT be"
        time.sleep(0.1)



@step("Verify the page is showing <number_of_items_text>")
def verify_showing_item_text(number_of_items_text):
    """Verifies the 'Showing' text at the bottom of tables has the correct text"""
    showing_items_text = APP.shared.get_showing_count_text()
    assert (
        f"{showing_items_text}" in showing_items_text
    ), f"The page is NOT showing expected number of items - {number_of_items_text}"

@step("Verify the table <table_name> is showing <number_of_items_text>")
def verify_table_showing_item_text(table_name, number_of_items_text):
    """Verifies the specified table's 'Showing' text is correct"""
    showing_items_text = APP.shared.get_table_showing_count_text(table_name)
    assert (
        f"{showing_items_text}" in showing_items_text
    ), f"The table '{table_name}' is NOT showing expected number of items - {number_of_items_text}"

@step("Verify the table <table_name> has a total of <number_of_items_text> items")
def verify_table_total_of_text(table_name, number_of_items_text):
    """Verifies the specified table's 'Total Of' text is correct"""
    total_of_items_text = APP.shared.get_table_item_count_text(table_name)
    assert (
        total_of_items_text == number_of_items_text
    ), f"The table '{table_name}' is NOT showing expected total of items - {number_of_items_text}"


@step("Verify the table header text is correct <table>")
def verify_table_header_text(table):
    """Verifies the table header has the correct text"""
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        table_header_text_by_column = APP.shared.get_table_header_text_by_column(v[1])
        # Remove new lines from input
        table_header_text_by_column = table_header_text_by_column.replace("\n", "")
        # Remove unwanted additional spaces between words from input
        table_header_text_by_column = re.sub(" +", " ", table_header_text_by_column)
        assert (
            f"{table_header_text_by_column}" == v[0]
        ), f"The table header column '{v[1]}' is showing text '{table_header_text_by_column}' when we expected text '{v[0]}'"

@step("Verify the table <table_name> header text is correct <table>")
def verify_table_header_text(table_name:str, table):
    """Verifies the table header has the correct text"""
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        table_header_text_by_column = APP.shared.get_table_header_text_by_column_in_specified_table(table_name, v[1])
        # Remove new lines from input
        table_header_text_by_column = table_header_text_by_column.replace("\n", "")
        # Remove unwanted additional spaces between words from input
        table_header_text_by_column = re.sub(" +", " ", table_header_text_by_column)
        assert (
            f"{table_header_text_by_column}" == v[0]
        ), f"In table '{table_name}': The table header column '{v[1]}' is showing text '{table_header_text_by_column}' when we expected text '{v[0]}'"

@step("Verify the vertical table <table_name> header text is correct <table>")
def verify_table_header_text(table_name:str, table):
    """Verifies the vertical table header has the correct text"""
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        table_header_text_by_row = APP.shared.get_table_header_text_by_row_in_specified_vertical_table(table_name, v[1])
        # Remove new lines from input
        table_header_text_by_row = table_header_text_by_row.replace("\n", "")
        # Remove unwanted additional spaces between words from input
        table_header_text_by_row = re.sub(" +", " ", table_header_text_by_row)
        assert (
            f"{table_header_text_by_row}" == v[0]
        ), f"In table '{table_name}': The table header column '{v[1]}' is showing text '{table_header_text_by_row}' when we expected text '{v[0]}'"

@step("Verify the table body text is correct <table>")
def verify_table_body_text(table):
    """Verifies the table body has the correct text"""
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        table_body_text_by_row_column = APP.shared.get_table_body_text_by_row_column(
            v[1], v[2]
        )
        # Remove new lines from input
        table_body_text_by_row_column = table_body_text_by_row_column.replace("\n", "")
        # Remove unwanted additional spaces between words from input
        table_body_text_by_row_column = re.sub(" +", " ", table_body_text_by_row_column)
        assert (
            f"{table_body_text_by_row_column}" == v[0]
        ), f"The table body row '{v[1]}' and column '{v[2]}' is showing text '{table_body_text_by_row_column}' when we expected text '{v[0]}'"

@step("Verify the table <table_name> body text is correct <table>")
def verify_table_body_text_in_specified_table(table_name:str, table):
    """
    In specified table, verifies the table body has the correct text.

    table_name - Table to verify
    v[0] - The expected text to appear
    v[1] - The row of the table
    v[2] - The column of the table
    """
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        """Verifies the table body has the correct text"""
        table_body_text_by_row_column = APP.shared.get_specified_table_body_text_by_row_column(
            table_name, v[1], v[2]
        )
        # Remove new lines from input
        table_body_text_by_row_column = table_body_text_by_row_column.replace("\n", "")
        # Remove unwanted additional spaces between words from input
        table_body_text_by_row_column = re.sub(" +", " ", table_body_text_by_row_column)
        assert (
            f"{table_body_text_by_row_column}" == v[0]
        ), f"In table '{table_name}': The table body row '{v[1]}' and column '{v[2]}' is showing text '{table_body_text_by_row_column}' when we expected text '{v[0]}'"


@step("Verify the table body tooltips are correct <table>")
def verify_table_body_tooltips_text(table):
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    """Verifies the table body has correct tooltips"""
    for k, v in enumerate(table):
        APP.shared.hover_table_body_by_row_column(v[1], v[2])
        is_tooltip_text_present = APP.shared.is_text_present(v[0])
        assert (
            is_tooltip_text_present
        ), f"Hovering over table body row '{v[1]}' and column '{v[2]}' does NOT produce the tooltip '{v[0]}' as we expect"

@step("Verify the table <table_name> body tooltips are correct <table>")
def verify_table_body_tooltips_text(table_name:str, table):
    """
    In specified table, verifies correct tooltip appears.

    table_name - Table to validate
    v[0] - The expected text to appear
    v[1] - The row of the table
    v[2] - The column of the table
    """

    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    """Verifies the table body has correct tooltips"""
    for k, v in enumerate(table):
        APP.shared.hover_table_body_by_row_column_in_specified_table(table_name, v[1], v[2])
        is_tooltip_text_present = APP.shared.is_text_present(v[0])
        assert (
            is_tooltip_text_present
        ), f"Hovering over table body row '{v[1]}' and column '{v[2]}' does NOT produce the tooltip '{v[0]}' as we expect"


@step("Verify the table <table_name> is displaying this information <table>")
def verify_table_is_displaying_text(table_name, table):
    """Verifies the table is displaying given text"""
    for k, v in enumerate(table):
        is_table_text_present = APP.shared.is_table_displaying_text(table_name, v[0])
        assert (
            is_table_text_present
        ), f"The table '{table_name}' is NOT displaying '{v[0]}'"

@step("Verify the table <table_name> is not displaying this information <table>")
def verify_table_is_not_displaying_text(table_name, table):
    """Verifies the table is not displaying given text"""
    for k, v in enumerate(table):
        is_table_text_present = APP.shared.is_table_displaying_text_no_wait(table_name, v[0])
        assert (
            is_table_text_present == False
        ), f"The table '{table_name}' is displaying '{v[0]}' when it should NOT"

@step("Verify the button <button_name> is disabled")
def verify_button_is_disabled(button_name: str):
    is_button_disabled = APP.shared.is_button_disabled(button_name)
    assert (
        is_button_disabled
    ), f"The button '{button_name}' is NOT disabled when it should be"


@step("Verify the button <button_name> is enabled")
def verify_button_is_disabled(button_name: str):
    is_button_disabled = APP.shared.is_button_disabled(button_name)
    assert (
        is_button_disabled == False
    ), f"The button '{button_name}' is disabled when it should NOT be"

@step("Verify the button with displayed text <button_name> is disabled")
def verify_named_button_is_disabled(button_name: str):
    is_button_disabled = APP.shared.is_named_button_disabled(button_name)
    assert (
        is_button_disabled
    ), f"The button '{button_name}' is disabled when it should NOT be"

@step("Verify the button <button_name> area is <expanded_or_collapsed>")
def verify_named_button_is_disabled(button_name: str, expanded_or_collapsed:str):
    is_button_expanded = APP.shared.is_button_area_expanded(button_name)
    if expanded_or_collapsed.lower() == "expanded":
        assert is_button_expanded == 'true', f"The button '{button_name}' is NOT expanded when it should be"
    elif expanded_or_collapsed.lower() == "collapsed":
        assert is_button_expanded == 'false', f"The button '{button_name}' is NOT collapsed when it should be"

@step("Wait for <data_testid> to be present on the page")
def wait_for_data_testid_to_be_visible_on_the_page(data_testid: str):
    """Waits for specified data-testid to be present on the page"""
    is_data_testid_visible = APP.shared.wait_for_data_testid_to_be_visible(data_testid)
    assert is_data_testid_visible, f"The data-testid '{data_testid}' is NOT present"


@step("Wait for loading spinner")
def wait_for_loading_spinner_generic_to_appear_then_disappear():
    """Waits for loading spinner to appear and disappear on the page"""
    # Wait for the spinner to appear.
    # If it does not, we move on; No need to fail the test.
    try:
        APP.shared.wait_for_loading_spinner_to_be_visible(15000)
        APP.shared.wait_for_loading_spinner_to_detatch()
    except:
        APP.shared.wait_for_loading_spinner_to_detatch()


@step("Wait for cohort bar case count loading spinner")
def wait_for_loading_spinner_cohort_bar_case_count_to_disappear():
    """Waits for cohort bar case count loading spinner to disappear on the page"""
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()


@step("Wait for table loading spinner")
def wait_for_loading_spinner_cohort_bar_case_count_to_disappear():
    """Waits for table loading spinner to disappear on the page"""
    APP.shared.wait_for_loading_spinner_table_to_detatch()

@step("Wait for loading spinners to disappear")
def wait_for_loading_spinners_to_disappear():
    """Waits for all loading spinners to detach"""
    APP.shared.wait_for_loading_spinners_to_detach()

@step("Wait for table body text to appear <table>")
def wait_for_table_body_text_to_appear(table):
    """Waits for specified table body text to appear"""
    # Wait for all possible loading spinners to detach before checking text
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        """
        v[0] - Text
        v[1] - Row
        v[2] - Column
        """
        APP.shared.wait_for_table_body_text_by_row_column(v[0], v[1], v[2])
        time.sleep(2)
        # Occasionally, the screen flickers where it shows the text we
        # are waiting for then it disappears for a moment. Checking for the
        # text twice should account for that.
        APP.shared.wait_for_table_body_text_by_row_column(v[0], v[1], v[2])
        time.sleep(2)
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()

@step("Wait for table <table_name> body text to appear <table>")
def wait_for_table_body_text_to_appear(table_name:str, table):
    """In specified table, waits for given table body text to appear"""
    # Wait for all possible loading spinners to detach before checking text
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    for k, v in enumerate(table):
        """
        v[0] - Text
        v[1] - Row
        v[2] - Column
        """
        APP.shared.wait_for_specified_table_body_text_by_row_column(table_name, v[0], v[1], v[2])
        time.sleep(2)
        # Occasionally, the screen flickers where it shows the text we
        # are waiting for then it disappears for a moment. Checking for the
        # text twice should account for that.
        APP.shared.wait_for_specified_table_body_text_by_row_column(table_name, v[0], v[1], v[2])
        time.sleep(2)
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()

@step("Is text <expected_text> present on the page")
def is_text_present_on_the_page(expected_text: str):
    """Verifies if expected text is on the page"""
    is_text_present = APP.shared.is_text_present(expected_text)
    assert is_text_present, f"The text '{expected_text}' is NOT present"


@step("Is text <expected_text> not present on the page")
def is_text_present_on_the_page(expected_text: str):
    """Verifies if text is no longer on the page as expected"""
    is_text_not_present = APP.shared.is_text_not_present(expected_text)
    assert (
        is_text_not_present
    ), f"The text '{expected_text}' is present when it should not"


@step("Is modal with text <expected_text> present on the page and <action>")
def is_modal_text_present_on_the_page(expected_text: str, action: str):
    """Waits for modal with specified text and optionally removes modal"""
    is_text_present = APP.shared.wait_for_text_in_temporary_message(
        expected_text, action
    )
    assert is_text_present, f"The text '{expected_text}' is NOT present in a modal"
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()


@step("Validate the message <message_id> displays the text <expected_text>")
def validate_message_id_text_is_present_on_the_page(
    message_id: str, expected_text: str
):
    """Verifies if specified data-testid message displays expected text"""
    is_text_present = APP.shared.is_message_id_text_present(message_id, expected_text)
    assert is_text_present, f"The text '{expected_text}' is NOT present"


@step("Collect these data portal statistics for comparison <table>")
def store_home_page_data_portal_statistics(table):
    """
    Stores data portal summary statistics for use in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'

    v[0] - The name of the home page statistic to collect
    v[1] - The name the statistic will be stored under
    """
    for k, v in enumerate(table):
        category_statistic = APP.home_page.get_data_portal_summary_statistic(v[0])
        data_store.spec[f"{v[1]}"] = category_statistic


@step("Collect button labels in table for comparison <table>")
def store_button_labels_in_tables_for_comparison(table):
    """
    Stores button label text for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'

    v[0] - The name of how the label will be stored
    v[1] - The row of the table
    v[2] - The column of the table
    """
    for k, v in enumerate(table):
        table_body_text_by_row_column = APP.shared.get_table_body_text_by_row_column(
            v[1], v[2]
        )
        data_store.spec[f"{v[0]}"] = table_body_text_by_row_column

@step("Collect button labels in table <table_name> for comparison <table>")
def store_button_labels_in_tables_for_comparison(table_name:str, table):
    """
    In specified table, stores button label text for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'

    table_name - Table to collect label from
    v[0] - The name of how the label will be stored
    v[1] - The row of the table
    v[2] - The column of the table
    """
    for k, v in enumerate(table):
        table_body_text_by_row_column = APP.shared.get_specified_table_body_text_by_row_column(
            table_name, v[1], v[2]
        )
        data_store.spec[f"{v[0]}"] = table_body_text_by_row_column

@step("Collect <text_id_to_collect> Item Count for comparison")
def store_text_item_count_for_comparison(text_id_to_collect:str):
    """
    Stores the given text data-testid item count for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'
    """
    data_store.spec[f"{text_id_to_collect} Item Count"] = APP.shared.get_text_by_data_testid(text_id_to_collect)

@step("Collect Cohort Bar Case Count for comparison")
def store_cohort_bar_case_count_for_comparison():
    """
    Stores current cohort bar case count for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'
    """
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    data_store.spec["Cohort Bar Case Count"] = APP.shared.get_cohort_bar_case_count()


@step("Collect <cohort_name> Case Count for comparison")
def store_cohort_bar_case_count_for_comparison(cohort_name: str):
    """
    Stores current cohort bar case count for comparison in future tests.
    Store the information using a key based on the cohort name.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'
    """
    data_store.spec[
        f"{cohort_name} Case Count"
    ] = APP.shared.get_cohort_bar_case_count()

@step("Collect table <table_name> Item Count for comparison")
def store_item_count_in_table_for_comparison(table_name:str):
    """
    Stores the total number of items in the given table for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'
    """
    data_store.spec[f"{table_name} Item Count"] = APP.shared.get_table_item_count_text(table_name)

@step("Collect case counts for the following filters for cohort <cohort_name> <table>")
def collect_case_counts_on_filters(cohort_name: str, table):
    """
    collect_case_counts_on_filters - Collect case count on filters.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'.
    :param cohort_name: Cohort Name we are collecting the information under
    :param v[0]: Filter Card Name
    :param v[1]: Filter we are collecting case count info on
    """
    for k, v in enumerate(table):
        # Expands list of filters to select if possible
        if APP.shared.is_show_more_or_show_less_button_visible_within_filter_card(
            v[0], "plus-icon"
        ):
            APP.shared.click_show_more_less_within_filter_card(
                v[0], "plus-icon"
            )
            time.sleep(0.5)

        case_count = (
            APP.shared.get_filter_selection_count(
                v[0], v[1]
            )
        )
        # Saves the case count under the filter, filter and cohort name
        data_store.spec[f"{v[0]}_{v[1]}_{cohort_name} Count"] = case_count

@step("The cohort bar case count should be <case_count>")
def is_cohort_bar_case_count_present_on_the_page(case_count: str):
    """Checks the cohort bar case count"""
    is_case_count_present = APP.shared.is_cohort_bar_case_count_present(case_count)
    assert (
        is_case_count_present
    ), f"The cohort bar is NOT displaying '{case_count}' cases"


@step("The cart should have <correct_file_count> files")
def is_cart_count_correct(correct_file_count: str):
    """Checks the cart file count in the upper-right corner of the data portal"""
    is_cart_count_correct = APP.shared.is_cart_count_correct(correct_file_count)
    assert (
        is_cart_count_correct
    ), f"The cart count is NOT displaying '{correct_file_count}'"


@step("Is data-testid button <data_testid> not present on the page")
def is_data_testid_not_present_on_the_page(data_testid: str):
    is_data_testid_present = APP.shared.is_data_testid_button_present(data_testid)
    assert (
        is_data_testid_present == False
    ), f"The data-testid '{data_testid}' IS present"


@step("Select <data_testid> on page")
def click_data_testid(data_testid: str):
    """Clicks specified data-testid"""
    APP.shared.click_data_testid(data_testid)


@step("Select <data_testid> a data-testid button")
def click_button_with_data_testid(data_testid: str):
    """Clicks specified data-testid button"""
    APP.shared.click_button_data_testid(data_testid)


@step("Select button <data_testid>")
def click_button_with_data_testid(data_testid: str):
    """Normalizes identifier, and clicks specified data-testid button"""
    data_testid = APP.shared.normalize_button_identifier(data_testid)
    APP.shared.click_button_data_testid(data_testid)


@step("Select <button_text_name>")
def click_button_with_displayed_text_name(button_text_name: str):
    """Selects a button based on displayed text"""
    APP.shared.click_button_with_displayed_text_name(button_text_name)

@step("Select the link <link_data_testid>")
def click_link_data_testid(link_data_testid: str):
    """Clicks a link with a data-testid"""
    APP.shared.click_link_data_testid(link_data_testid)

@step("Select the following checkboxes <table>")
def click_checkboxes(table):
    """Selects checkboxes given from table exactly as typed"""
    for k, v in enumerate(table):
        APP.shared.click_checkbox(v[0])
        time.sleep(0.1)

@step("Select checkbox <checkbox_name>")
def click_checkbox_normalize_ident(checkbox_name:str):
    """Normalizes the checkbox name, then attempts to check"""
    APP.shared.click_checkbox_normalize_ident(checkbox_name)
    time.sleep(0.1)

@step("Select the radio button <button_name>")
def click_radio_data_testid_button(button_name:str):
    APP.shared.click_radio_data_testid_button(button_name)

@step("Select the following radio buttons <table>")
def click_radio_buttons(table):
    for k, v in enumerate(table):
        APP.shared.click_radio_button(v[0])
        time.sleep(0.1)

@step("Select create or save in cohort modal")
def click_create_or_save_in_cohort_modal():
    """Clicks 'Create' or 'Save' in cohort modal"""
    APP.shared.click_create_or_save_button_in_cohort_modal()


@step("Select or deselect these options from the table column selector <table>")
def click_create_or_save_in_cohort_modal(table):
    """
    Clicks table column selector button.
    In the column selector pop-up modal that appears, it clicks the specified switch.
    """
    APP.shared.click_column_selector_button()
    for k, v in enumerate(table):
        APP.shared.click_switch_for_column_selector(v[0])
    APP.shared.click_column_selector_button()

@step("In table <table_name> select or deselect these options from the table column selector <table>")
def click_create_or_save_in_cohort_modal(table_name:str, table):
    """
    In specified table, clicks table column selector button.
    In the column selector pop-up modal that appears, it clicks the specified switch.
    """
    APP.shared.click_column_selector_button_in_specified_table(table_name)
    for k, v in enumerate(table):
        APP.shared.click_switch_for_column_selector(v[0])
    APP.shared.click_column_selector_button_in_specified_table(table_name)

@step("In table <table_name> restore default column selector options")
def click_restore_column_selector_options(table_name:str):
    """In specified table, clicks reset column selector button"""
    APP.shared.click_column_selector_button_in_specified_table(table_name)
    APP.shared.click_reset_column_select_options()
    APP.shared.click_column_selector_button_in_specified_table(table_name)

@step("Select <text> from dropdown menu")
def click_dropdown_text_option(text:str):
    """Selects text option from dropdown menu with data-testid dropdown-menu-options"""
    APP.shared.click_text_option_from_dropdown_menu(text)

@step("Select partial text <text> from dropdown menu")
def click_dropdown_text_option(text:str):
    """Selects text option from dropdown menu with data-testid dropdown-menu-options"""
    APP.shared.click_has_text_option_from_dropdown_menu(text)

@step("Change number of entries shown in the table to <number_of_entries>")
def change_number_of_entries_shown(change_number_of_entries_shown: str):
    """
    Changes number of entries shown in the table using the show entries button,
    and selecting an option from the dropdown list.
    """
    APP.shared.change_number_of_entries_shown(change_number_of_entries_shown)
    APP.shared.wait_for_loading_spinners_to_detach()

@step("Change number of entries shown in the table <table_name> to <number_of_entries>")
def change_number_of_entries_shown_in_specified_table(table_name:str, change_number_of_entries_shown: str):
    """
    In specified table, changes number of entries shown in the table using the show entries button,
    and selecting an option from the dropdown list.
    """
    APP.shared.change_number_of_entries_shown_in_specified_table(table_name, change_number_of_entries_shown)
    APP.shared.wait_for_loading_spinners_to_detach()

@step("Perform action and validate modal text <table>")
def click_named_button_in_modal_and_wait_for_temp_message_text(table):
    """
    click_named_button_wait_for_message_text clicks a button by its displayed name in a modal,
    validates text in a temporary pop-up modal message, and either clicks the 'x' to remove the temp message,
    or does nothing to let the message persist.

    :param v[0]: The name of the button to be clicked
    :param v[1]: The text in the temporary message that we are waiting for
    :param v[2]: Input of "Removal Modal" will remove the temp message, otherwise we let it persist
    """
    APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
    for k, v in enumerate(table):
        APP.shared.click_button_in_modal_with_displayed_text_name(v[0])
        is_cohort_message_present = APP.cohort_bar.wait_for_text_in_temporary_message(
            v[1], v[2]
        )
        assert is_cohort_message_present, f"The text '{v[1]}' is NOT present"
        # Need to let the page load after our actions here.
        # Automation moves too quickly in the cohort bar section.
        time.sleep(1)
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()

@step("Remove message if one is present")
def remove_message_if_one_is_present():
    if APP.shared.is_temporary_message_present():
        APP.shared.click_close_temporary_message()

@step("Clear active cohort filters")
def clear_active_cohort_filters():
    # Clicks the 'clear all' button in the cohort query area
    APP.shared.clear_active_cohort_filters()


@step("Undo Action")
def click_undo_in_message():
    """Clicks 'undo' in a modal message"""
    APP.shared.click_undo_in_message()


@step("Set this as your current cohort")
def click_undo_in_message():
    """Clicks 'Set this as your current cohort' in a modal message"""
    APP.shared.click_set_as_current_cohort_in_message()
    APP.shared.wait_for_loading_spinners_to_detach()


# These 3 functions are for filter cards (like on projects page).
# The filter cards depend on a specific data-testid "filters-facets" that
# is common across a couple of filter sets in the data portal.
@step("Make the following selections on a filter card <table>")
def filter_card_selections(table):
    """Trio of actions for the filter cards and filters on the projects page"""
    for k, v in enumerate(table):
        # Expands list of filters to select if possible
        if APP.shared.is_show_more_or_show_less_button_visible_within_filter_card(
            v[0], "plus-icon"
        ):
            APP.shared.click_show_more_less_within_filter_card(v[0], "plus-icon")
        APP.shared.make_selection_within_filter_group(v[0], v[1])
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()
        time.sleep(0.1)


@step("Perform the following actions on a filter card <table>")
def perform_filter_card_action(table):
    for k, v in enumerate(table):
        APP.shared.perform_action_within_filter_card(v[0], v[1])
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()
        time.sleep(0.1)

@step("Flip the switch on filter card <filter_card_name>")
def flip_switch_filter_card(filter_card_name:str):
    APP.shared.flip_switch_in_filter_card(filter_card_name)
    time.sleep(1)
    APP.shared.wait_for_loading_spinners_to_detach()

@step("Search in a filter card <table>")
def search_in_filter_card(table):
    for k, v in enumerate(table):
        APP.shared.type_in_filter_card_search_text_area(v[0], v[1], v[2])
        time.sleep(0.1)

@step("Expand or contract a filter <table>")
def click_show_more_or_show_less(table):
    for k, v in enumerate(table):
        APP.shared.click_show_more_less_within_filter_card(v[0], v[1])

@step("Expand dropdown in graph <graph_id>")
def click_dropdown_button_in_graph(graph_id):
    APP.shared.click_button_dropdown_in_graph(graph_id)
    time.sleep(0.5)

@step("Select value from table by row and column <table>")
def select_table_value_by_row_column(table):
    """
    Selects values from tables by giving a row and column
    Row and Column indexing begins at '1'
    """
    for k, v in enumerate(table):
        APP.shared.select_table_by_row_column(v[0], v[1])
        time.sleep(1)
        # In Mutation Frequency, selecting items in the table can take a
        # long time to load. They can also load and spin at different times
        # in different places (e.g the cohort case count, table, graphs, etc.)
        # So we have an abundance of waits.
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()

@step("Select value from table <table_name> by row and column <table>")
def select_table_value_by_row_column(table_name:str, table):
    """
    In specified table, selects values from tables by giving a row and column
    Row and Column indexing begins at '1'
    """
    for k, v in enumerate(table):
        APP.shared.select_specified_table_by_row_column(table_name, v[0], v[1])
        time.sleep(1)
        # In Mutation Frequency, selecting items in the table can take a
        # long time to load. They can also load and spin at different times
        # in different places (e.g the cohort case count, table, graphs, etc.)
        # So we have an abundance of waits.
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()

@step("Enter text <text> in the <aria_label> search bar")
def send_text_into_search_bar(text: str, aria_label: str):
    """Sends text into search bar based on its aria_label"""
    APP.shared.send_text_into_search_bar(text, aria_label)


@step("Enter <text> in the text box <text_box_name>")
def send_text_into_text_box(text: str, text_box_name: str):
    """Sends text into a data-testid text box"""
    APP.shared.send_text_into_text_box(text, text_box_name)


@step("In table <table_name>, search the table for <text>")
def send_text_into_specified_table_search_bar(table_name: str, text: str):
    """Sends text into specified table search bar"""
    APP.shared.send_text_into_specified_table_search_bar(table_name, text)
    time.sleep(1)
    # In Mutation Frequency, searching in the table can take a
    # long time to load. They can also load and spin at different times
    # in different places (e.g the cohort case count, table, graphs, etc.)
    # So we have an abundance of waits.
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.keyboard_press("Enter")
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()


@step("Search the table for <text>")
def send_text_into_table_search_bar(text: str):
    """Sends text into a table search bar"""
    APP.shared.send_text_into_table_search_bar(text)
    time.sleep(1)
    # In Mutation Frequency, searching in the table can take a
    # long time to load. They can also load and spin at different times
    # in different places (e.g the cohort case count, table, graphs, etc.)
    # So we have an abundance of waits.
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.keyboard_press("Enter")
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.shared.wait_for_loading_spinner_to_detatch()


@step("Quick search for <text> and go to its page")
def quick_search_and_click(text: str):
    """
    Sends text into the quick search bar in the upper-right corner of the data portal.
    Then clicks the result in the search result area. Best used with a UUID.
    """
    APP.shared.quick_search_and_click(text)


@step("Quick search for <text>")
def global_quick_search(text: str):
    """
    Sends text into the quick search bar in the upper-right corner of the data portal.
    """
    APP.shared.global_quick_search(text)


@step(
    "Validate the quick search bar result in position <result_in_list> of the result list has the text <text>"
)
def validate_global_quick_search_result_text(result_in_list: str, text: str):
    """
    Specifies a result from the quick search bar result list. Validates expected text is present.
    """
    APP.shared.validate_global_quick_search_result_text(result_in_list, text)


@step("Select the quick search bar result in position <result_in_list>")
def quick_search(result_in_list: str):
    """
    Specifies a result from the quick search bar result list. Clicks that result.
    """
    APP.shared.click_global_quick_search_result(result_in_list)


@step("Name the cohort <cohort_name>")
def name_cohort(cohort_name: str):
    APP.shared.send_text_into_text_box(cohort_name, "Name Input Field")


@step(
    "These links on the <page_name> should take the user to correct page in a new tab <table>"
)
def click_nav_item_check_text_in_new_tab(page_name: str, table):
    """
    Performs an action to open a new tab.
    Then, checks for expected text on the new tab to indicate it opened correctly.
    """

    # Banners can interfere with this test. Dismiss any banner before beginning.
    if APP.shared.is_text_present("Dismiss"):
        APP.shared.click_button_with_displayed_text_name("Dismiss")
    for k, v in enumerate(table):
        new_tab = APP.shared.perform_action_handle_new_tab(page_name, v[0])
        is_text_visible = APP.shared.is_text_visible_on_new_tab(new_tab, v[1])
        assert (
            is_text_visible
        ), f"After click on '{v[0]}', the expected text '{v[1]}' in NOT present"
        new_tab.close()

@step(
    "In table <table_name> these selections should take the user to correct page in a new tab <table>"
)
def click_table_by_row_column_check_text_in_new_tab(table_name: str, table):
    """
    click_table_by_row_column_check_text_in_new_tab clicks in a table to open a new tab,
    and validates the text present in the new tab.

    :param table_name: Specifies what table to click on.
    :param v[0]: Row number to click.
    :param v[1]: Column number to click.
    :param v[2]: Text to check.
    """
    for k, v in enumerate(table):
        new_tab = APP.shared.click_in_table_handle_new_tab(table_name, v[0], v[1])
        is_text_visible = APP.shared.is_text_visible_on_new_tab(new_tab, v[2])
        assert (
            is_text_visible
        ), f"After click in table '{table_name}', row: '{v[0]}' and column: '{v[1]}' the expected text '{v[2]}' in NOT present"
        new_tab.close()

@step(
    "Check that <var_to_check> cookie is accessible using Javascript and that it's generated using uuid version <ver>"
)
def check_if_cookie_accessible(var_to_check: str, ver: int):
    cookie = APP.driver.evaluate("()=>document.cookie")
    var_to_check += "="

    # check to see if gdc_context_id exists in cookie
    assert var_to_check in cookie

    start_value = cookie[cookie.index(var_to_check) + len(var_to_check) :]
    gdc_context_id = start_value[: start_value.find(";")]
    start_value = cookie[cookie.index(var_to_check) + len(var_to_check) :]
    gdc_context_id = start_value[: start_value.find(";")]

    # check if the gdc_context_id is version 4
    assert UUID(gdc_context_id).version == int(ver)
