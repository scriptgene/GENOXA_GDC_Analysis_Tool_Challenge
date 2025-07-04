import time

from getgauge.python import step, before_spec

from ..app import GDCDataPortalV2App
from ....base.webdriver import WebDriver

from getgauge.python import data_store


@before_spec
def start_app():
    global APP
    APP = GDCDataPortalV2App(WebDriver.page)


@step("Select <button_name> on the Repository page")
def select_repository_page_button(button_name: str):
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    APP.repository_page.click_repository_page_button(button_name)


@step("Select <button_name> on the Image Viewer page")
def select_image_viewer_page_button(button_name: str):
    APP.repository_page.click_image_viewer_page_data_testid(button_name)
    APP.repository_page.wait_for_loading_spinner_to_detatch()


@step("Select case or slide <data_testid> on the Image Viewer page")
def select_case_slide_on_image_viewer(data_testid: str):
    APP.repository_page.click_image_viewer_page_case_or_slide(data_testid)
    APP.repository_page.wait_for_loading_spinner_to_detatch()


@step("Select file filter, <filter_name>, nth: <nth>")
def select_file_filter_and_validate(filter_name: str, nth: int):
    repository = APP.repository_page
    try:
        repository.click_button(filter_name)
    except:
        repository.select_nth_file_filters_result(int(nth) - 1)

@step("Add the following files to the cart on the Repository page <table>")
def add_files_to_cart(table):
    """
    add_files_to_cart Searches the repository table for a File UUID and adds
    that file to the cart.
    :param v[0]: File UUID to add
    """
    APP.shared.click_column_selector_button_in_specified_table("repository")
    APP.shared.click_switch_for_column_selector("file_uuid")
    APP.shared.click_column_selector_button_in_specified_table("repository")
    for k, v in enumerate(table):
        APP.shared.send_text_into_specified_table_search_bar("repository",v[0])
        APP.shared.wait_for_specified_table_body_text_by_row_column("repository", v[0], 1, 2)
        APP.shared.select_specified_table_by_row_column("repository", 1, 1)

@step("Collect file counts for the following filters on the Repository page <table>")
def collect_file_counts_on_filters(table):
    """
    collect_file_counts_on_filters Collect file count on filters on the repository page.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'.
    :param v[0]: Facet Card Name
    :param v[1]: Filter we are collecting file count info on
    """
    for k, v in enumerate(table):
        # Expands list of filters to select if possible
        if APP.repository_page.is_show_more_or_show_less_button_visible_within_filter_card_repository(
            v[0], "plus-icon"
        ):
            APP.repository_page.click_show_more_less_within_filter_card_repository(
                v[0], "plus-icon"
            )

        file_count = (
            APP.repository_page.get_file_count_from_filter_within_facet_group(
                v[0], v[1]
            )
        )
        # Saves the file count under the facet and filter
        data_store.spec[f"{v[0]}_{v[1]} Repository Count"] = file_count

@step("Collect <item> Count on the Repository page")
def store_count_repository_for_comparison(item: str):
    """
    Stores specified item count for comparison in future tests.
    Pairs with the test 'verify_compared_statistics_are_equal_or_not_equal'

    :param item: The item to collect the count of. Options: Files, Cases, or Size
    """

    # Everything on the repository page could load, except for the file, case, size count.
    # So, we check to see if the text displays '--' which means it's still loading.
    # We check for the text to display an actual number for approximately 15 seconds.
    # If no number appears, we return the '--' as is, and the test will fail as intended.
    item_count = APP.repository_page.get_repository_table_item_count(item)
    retry_counter = 0
    while item_count == "--":
        time.sleep(1)
        item_count = APP.repository_page.get_repository_table_item_count(item)
        retry_counter = retry_counter+1
        if retry_counter >= 15:
            break
    data_store.spec[f"{item} Count Repository Page"] = item_count

@step("Verify filters have correct access level <table>")
def verify_filters_access_level(table):
    """
    verify_filters_access_level Selects a filter and then compares its file count
    to the file count of given access level. To be used when a filter should have
    all Open Access or all Controlled Access.

    :param v[0]: Facet card
    :param v[1]: Filter to select and compare
    :param v[2]: Access Level
    """
    for k, v in enumerate(table):
        if APP.repository_page.is_show_more_or_show_less_button_visible_within_filter_card_repository(
            v[0], "plus-icon"
        ):
            APP.repository_page.click_show_more_less_within_filter_card_repository(
                v[0], "plus-icon"
            )

        APP.shared.make_selection_within_filter_group(v[0], v[1])
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()
        time.sleep(0.1)

        filter_case_count = (
            APP.repository_page.get_file_count_from_filter_within_facet_group(
                v[0], v[1]
            )
        )

        access_case_count = (
            APP.repository_page.get_file_count_from_filter_within_facet_group(
                "Access", v[2]
            )
        )

        assert (
            filter_case_count == access_case_count
        ), f"The filter '{v[1]}' has value '{filter_case_count}' and access level '{v[2]}' has value '{access_case_count}' that does NOT match"


        APP.shared.perform_action_within_filter_card(v[0], "clear selection")
        APP.shared.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        APP.shared.wait_for_loading_spinner_table_to_detatch()
        APP.shared.wait_for_loading_spinner_to_detatch()
        time.sleep(0.1)

@step("Verify that the following default filters are displayed in order <table>")
def default_filters(table):
    repository = APP.repository_page
    actual_filters_order = repository.get_filter_facet_names()
    expected_default_filters_in_order = []
    for row, value in enumerate(table):
        expected_default_filters_in_order += value
    assert (
        len(actual_filters_order) == len(expected_default_filters_in_order)
        and actual_filters_order == expected_default_filters_in_order
    ), (
        f"Default filters sizes do not match with expected or they are"
        f" not displayed in order.\nActual: {actual_filters_order}"
        f"\nExpected: {expected_default_filters_in_order}"
    )

@step("Verify <count> items on Add a Custom Filter filter list")
def verify_file_filter_list_count(count: int):
    actual_count = APP.repository_page.get_file_filter_list_count()
    assert (
        int(count) == actual_count
    ), f"Add a Custom Filters list count mismatch.\nExpected: {count}\nActual: {actual_count}"


@step("Verify file filter names do not start with <pattern>")
def verify_file_filter_names_are_appropriate(pattern):
    file_filter_names = APP.repository_page.get_file_filter_names()
    fails = []
    for name in file_filter_names:
        if name.startswith(pattern):
            fails.append({name: name.startswith(pattern)})
    assert (
        not fails
    ), f"Some files are starting with `{pattern}`!\nFilter names: {fails}"


@step("Verify that the file filter, <filter_name>, has been applied")
def verify_file_filter_applied(filter_name: str):
    repository = APP.repository_page
    expected_filter_name = repository.get_custom_filter_facet_as_applied(filter_name)
    actual_filter_name = repository.get_filter_facet_names()[0]
    assert (
        expected_filter_name == actual_filter_name
    ), f"Custom filter not found in facets.\nExpected: {expected_filter_name}\nActual: {actual_filter_name}"


@step("Verify the slide image is visible")
def verify_slide_image_is_visible():
    is_slide_image_visible = APP.repository_page.is_slide_image_visible()
    assert is_slide_image_visible, f"The slide image is NOT visible"


@step("Verify the slide image viewer is showing <number_of_cases> cases")
def verify_slide_image_viewer_case_count(number_of_cases):
    APP.shared.wait_for_loading_spinner_table_to_detatch()
    showing_cases_slide_image_text = (
        APP.repository_page.get_image_viewer_showing_cases_text()
    )
    assert (
        f"Showing {number_of_cases}" in showing_cases_slide_image_text
    ), f"The slide image viewer page is not showing expected number of cases - {number_of_cases}"


@step("Verify details fields and values <table>")
def verify_details_fields_and_values(table):
    """On a slide image, details pop-up, verifies that given fields and values are present"""
    for k, v in enumerate(table):
        field_present = APP.repository_page.is_detail_field_present(v[0])
        assert (
            field_present
        ), f"Expected field '{v[0]}' is NOT present in details section"
        value_present = APP.repository_page.is_detail_value_present(v[0], v[1])
        assert (
            value_present
        ), f"Expected value '{v[1]}' for field '{v[0]}' is NOT present in details section"


@step("Search for file filter, <filter_name>")
def search_for_filter_and_searchbox_content(filter_name: str):
    APP.repository_page.search_file_filters(filter_name)
    search_box_content = APP.repository_page.get_search_box_entry()
    assert (
        filter_name == search_box_content
    ), f"Expected value not found in search box entry.\nExpected: {filter_name}\nActual: {search_box_content}"


@step("Search for <image_viewer_search> on the Image Viewer page")
def search_image_viewer(image_viewer_search: str):
    APP.repository_page.search_image_viewer(image_viewer_search)
    APP.repository_page.wait_for_loading_spinner_to_detatch()


@step("Remove search filter <search_filter> on the Image Viewer page")
def remove_slide_image_viewer_search_filter(search_filter: str):
    APP.repository_page.remove_slide_image_viewer_search_filter(search_filter)
    APP.repository_page.wait_for_loading_spinner_to_detatch()
