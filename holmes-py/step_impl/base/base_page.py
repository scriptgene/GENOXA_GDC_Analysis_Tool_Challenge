import time
import re

from typing import List
from step_impl.base.webdriver import WebDriver


class GenericLocators:
    TEXT_IDENT = lambda text: f"text={text} >> nth=0"
    TEXT_IN_PARAGRAPH = lambda text: f'p:has-text("{text}") >> nth=0'

    UNDO_BUTTON_IN_TEMP_MESSAGE = 'span:text("Undo")'
    SET_AS_CURRENT_COHORT_IN_TEMP_MESSAGE = (
        'span:text("Set this as your current cohort.")'
    )

    BUTTON_CLOSE_NOTIFICATION = '[aria-label="Close notification"] >> nth=0'
    BUTTON_CLOSE_MODAL = 'button[aria-label="Close Modal"]'

    LOADING_SPINNER_GENERIC = '[data-testid="loading-spinner"] >> nth=0'
    LOADING_SPINNER_COHORT_BAR_CASE_COUNT = (
        '[data-testid="loading-spinner-cohort-case-count"] >> nth=0'
    )
    LOADING_SPINNER_TABLE = '[data-testid="loading-spinner-table"] >> nth=0'

    # This locator allows you to send in a specific case count and check if it is there
    COHORT_BAR_CASE_COUNT = (
        lambda case_count: f'[data-testid="expandcollapseButton"] >> text="{case_count}"'
    )
    # This locator allows you to grab the case count text for further testing
    TEXT_COHORT_BAR_CASE_COUNT = (
        f'[data-testid="expandcollapseButton"] >> [data-testid="button-cases-cohort-bar"] >> span:nth-child(1)'
    )
    CART_IDENT = '[data-testid="button-header-cart"]'

    BUTTON_CLEAR_ACTIVE_COHORT_FILTERS = (
        '[data-testid="button-clear-all-cohort-filters"]'
    )
    TEXT_NO_ACTIVE_COHORT_FILTERS = '[data-testid="text-no-active-cohort-filter"]'

    BUTTON_CREATE_OR_SAVE_COHORT_IN_MODAL = '[data-testid="button-save-name"]'

    TEXT_BOX_IDENT = lambda text_box: f'[data-testid="textbox-{text_box}"]'
    SEARCH_BAR_ARIA_IDENT = lambda aria_label: f'[aria-label="{aria_label}"]'
    SEARCH_BAR_TABLE_IDENT = '[data-testid="textbox-table-search-bar"] >> nth=0'
    SEARCH_BAR_IN_SPECIFIED_TABLE_IDENT = (
        lambda table_specified: f'[data-testid="table-{table_specified}"] >> [data-testid="textbox-table-search-bar"] >> nth=0'
    )
    QUICK_SEARCH_BAR_IDENT = '[data-testid="textbox-quick-search-bar"]'
    QUICK_SEARCH_BAR_FIRST_RESULT = '[data-testid="text-search-result"] >> nth=0'
    QUICK_SEARCH_BAR_NUMBERED_RESULT = (
        lambda result_in_list: f'[data-testid="text-search-result"] >> nth={result_in_list}'
    )
    QUICK_SEARCH_BAR_RESULT_TEXT = (
        lambda result_in_list, text: f'[data-testid="text-search-result"] >> nth={result_in_list} >> text={text}'
    )

    RADIO_BUTTON_IDENT = lambda radio_name: f'//input[@id="{radio_name}"]'
    RADIO_DATA_TESTID_BUTTON_IDENT = lambda radio_name: f'[data-testid="radio-{radio_name}"]'
    CHECKBOX_IDENT = (
        lambda checkbox_id: f'[data-testid="checkbox-{checkbox_id}"]'
    )

    DATA_TEST_ID_IDENT = lambda id: f'[data-testid="{id}"]'
    DATA_TESTID_BUTTON_IDENT = (
        lambda data_testid: f'[data-testid="button-{data_testid}"]'
    )
    DATA_TESTID_TEXT_IDENT = lambda data_testid: f'[data-testid="text-{data_testid}"]'
    DATA_TESTID_LINK_IDENT = lambda data_testid: f'[data-testid="link-{data_testid}"]'

    BUTTON_BY_DISPLAYED_TEXT = (
        lambda button_text_name: f'button:has-text("{button_text_name}") >> nth=0'
    )
    BUTTON_IN_MODAL_BY_DISPLAYED_TEXT = (
        lambda button_text_name: f'section[role="dialog"] >> button:has-text("{button_text_name}") >> nth=0'
    )
    BUTTON_A_BY_TEXT_IDENT = (
        lambda button_text_name: f'a:has-text("{button_text_name}") >> nth=0'
    )
    BUTTON_IN_FOOTER_BY_TEXT_IDENT = (
        lambda button_text_name: f'footer >> a:has-text("{button_text_name}") >> nth=0'
    )

    TABLE_ITEM_COUNT_IDENT = (
        lambda table_name: f'[data-testid="table-{table_name}"] >> [data-testid="text-total-item-count"]'
    )
    TABLE_TEXT_IDENT = (
        lambda table_name, table_text: f'[data-testid="table-{table_name}"] >> text={table_text} >> nth=0'
    )
    TABLE_AREA_TO_SELECT = (
        lambda row, column: f"tr:nth-child({row}) > td:nth-child({column}) > * >> nth=0"
    )
    TABLE_AREA_TO_COLLECT_IN_SPECIFIED_TABLE = (
        lambda table_name, row, column: f"[data-testid='table-{table_name}'] >> tr:nth-child({row}) > td:nth-child({column}) >> nth=0"
    )
    TABLE_AREA_TO_SELECT_IN_SPECIFIED_TABLE = (
        lambda table_name, row, column: f"[data-testid='table-{table_name}'] >> tr:nth-child({row}) > td:nth-child({column}) > * >> nth=0"
    )
    TABLE_AREA_TO_CLICK = (
        lambda row, column: f"tr:nth-child({row}) > td:nth-child({column}) > * > * >> nth=0"
    )
    TABLE_AREA_TO_CLICK_IN_SPECIFIED_TABLE = (
        lambda table_name, row, column: f"[data-testid='table-{table_name}'] >> tr:nth-child({row}) > td:nth-child({column}) > * > * >> nth=0"
    )
    TABLE_CHECKBOX_TO_CLICK_IN_SPECIFIED_TABLE = (
        lambda table_name, row, column: f"[data-testid='table-{table_name}'] >> tr:nth-child({row}) > td:nth-child({column}) >> [type='checkbox']"
    )
    TABLE_TEXT_TO_WAIT_FOR = (
        lambda text, row, column: f'tr:nth-child({row}) > td:nth-child({column}) > * >> nth=0 >> text="{text}"'
    )
    TABLE_TEXT_TO_WAIT_FOR_IN_SPECIFIED_TABLE = (
        lambda table_name, text, row, column: f'[data-testid="table-{table_name}"] >> tr:nth-child({row}) > td:nth-child({column}) > * >> nth=0 >> text="{text}"'
    )
    TEXT_SPECIFIED_TABLE_HEADER = lambda table_name, column: f'[data-testid="table-{table_name}"] >> tr > th:nth-child({column}) >> nth=0'
    TEXT_SPECIFIED_VERTICAL_TABLE_HEADER = lambda table_name, row: f'[data-testid="table-{table_name}"] >> tr:nth-child({row}) > th >> nth=0'
    TEXT_TABLE_HEADER = lambda column: f"tr > th:nth-child({column}) >> nth=0"
    TEXT_DROPDOWN_MENU_OPTION = (
        lambda dropdown_option: f'[data-testid="dropdown-menu-options"] >> text="{dropdown_option}"'
    )
    TEXT_HAS_DROPDOWN_MENU_OPTION = (
        lambda dropdown_option: f'[data-testid="dropdown-menu-options"] >> button:has-text("{dropdown_option}")'
    )
    BUTTON_TEXT_DROPDOWN_MENU_OPTION = (
        lambda dropdown_option: f'[data-testid="dropdown-menu-options"] >> text="{dropdown_option}" >> ..'
    )
    BUTTON_TEXT_DROPDOWN_IMAGE_DATA = (
        lambda dropdown_option: f'[data-testid="list-download-image-or-data-dropdown"] >> text="{dropdown_option}"'
    )

    BUTTON_COLUMN_SELECTOR = '[data-testid="button-column-selector-box"]'
    BUTTON_COLUMN_SELECTOR_IN_TABLE = lambda table_name: f'[data-testid="table-{table_name}"] >> [data-testid="button-column-selector-box"] >> nth=0'
    SWITCH_COLUMN_SELECTOR = (
        lambda switch_name: f'[data-testid="column-selector-popover-modal"] >> [data-testid="column-selector-row-{switch_name}"] label div >> nth=0'
    )
    BUTTON_RESET_COLUMN_SELECTIONS = '[data-testid="restore-default-icon"]'

    BUTTON_CUSTOM_FILTER_IN_PROPERTIES_TABLE = lambda button_name: f'[data-testid="list-file-filters"] >> text="{button_name}"'
    FILTER_GROUP_IDENT = (
        lambda group_name: f'[data-testid="filters-facets"] >> div:text-is("{group_name}")'
    )
    FILTER_GROUP_CUSTOM_FILTER_TEXT_IDENT = (
        lambda group_name, filter_text: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> text="{filter_text}"'
    )
    FILTER_GROUP_SELECTION_IDENT = (
        lambda group_name, selection: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> [data-testid="checkbox-{selection}"]'
    )
    FILTER_GROUP_SELECTION_COUNT_IDENT = (
        lambda group_name, selection: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> [data-testid="text-{selection}"]'
    )
    FILTER_GROUP_ACTION_IDENT = (
        lambda group_name, action: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> button[aria-label="{action}"]'
    )
    FILTER_GROUP_FLIP_SWITCH_IDENT = (
        lambda group_name: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> [data-testid="toggle-area"] >> div >> nth=0'
    )
    FILTER_GROUP_SHOW_MORE_LESS_IDENT = (
        lambda group_name, more_or_less: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> button[data-testid="{more_or_less}"]'
    )
    FILTER_GROUP_TEXT_AREA_IDENT = (
        lambda group_name, textbox_id: f'[data-testid="filters-facets"] >> div:has-text("{group_name}") >> [data-testid="textbox-{textbox_id}"]'
    )

    SHOWING_NUMBER_OF_ITEMS_IN_TABLE = lambda table_specified: f'[data-testid="table-{table_specified}"] >> [data-testid="text-showing-count"]'
    SHOWING_NUMBER_OF_ITEMS = '[data-testid="text-showing-count"]'

    BUTTON_ENTRIES_SHOWN = '[data-testid="button-show-entries"]'
    BUTTON_ENTRIES_SHOWN_IN_SPECIFIED_TABLE = lambda table_name: f'[data-testid="table-{table_name}"] >> [data-testid="button-show-entries"]'

    BUTTON_IN_GRAPH = lambda graph_name: f'[data-testid="graph-{graph_name}"] >> [data-testid="button-download-image-or-data"]'
    BUTTON_IN_MOST_FREQUENT_SOMATIC_MUTATIONS_TABLE = lambda button_name: f'[data-testid="table-most-frequent-somatic-mutations"] >> text="{button_name}"'

    DROPDOWN_LIST_CHANGE_NUMBER_OF_ENTRIES_SHOWN = (
        lambda number_of_entries: f'[data-testid="area-show-number-of-entries"] >> text="{number_of_entries}"'
    )
    DROPDOWN_LIST_CHANGE_NUMBER_OF_ENTRIES_SHOWN_IN_SPECIFIED_TABLE = (
        lambda table_name, number_of_entries: f'[data-testid="table-{table_name}"] >> [data-testid="area-show-number-of-entries"] >> text="{number_of_entries}"'
    )


class BasePage:
    def __init__(self, driver) -> None:
        self.driver = driver

    def goto(self, url):
        self.driver.goto(url)

    # Force: Whether to bypass the actionability checks
    def click(self, locator, force=False, timeout=45000):
        self.wait_until_locator_is_visible(locator)
        self.driver.locator(locator).click(force=force, timeout=timeout)

    def hover(self, locator, force=False):
        """Hover over given locator"""
        self.driver.locator(locator).hover(force=force)

    def get_text(self, locator):
        return self.driver.locator(locator).text_content()

    def get_inner_text(self, locator):
        """Returns only what a user can see on the page"""
        return self.driver.locator(locator).inner_text()

    def get_input_value(self, locator):
        return self.driver.locator(locator).input_value()

    def get_attribute(self, locator, name: str):
        return self.driver.locator(locator).get_attribute(name)

    def get_count(self, locator):
        return self.driver.locator(locator).count()

    def is_checked(self, locator):
        return self.driver.locator(locator).is_checked()

    def is_visible(self, locator):
        return self.driver.locator(locator).is_visible()

    def is_disabled(self, locator):
        """Returns if the locator has the attribute 'disabled'"""
        return self.driver.locator(locator).is_disabled()

    def is_enabled(self, locator):
        """Returns if the locator is enabled"""
        return self.driver.locator(locator).is_enabled()

    def send_keys(self, locator, text):
        return self.driver.locator(locator).fill(text)

    def keyboard_press(self, action):
        return self.driver.keyboard.press(action)

    def scroll_into_view_if_needed(self, locator):
        self.driver.locator(locator).scroll_into_view_if_needed()

    def normalize_button_identifier(self, button_name: str) -> str:
        """Takes BDD spec file input and converts it to the ID formatting in the data portal"""
        return button_name.lower().replace(" ", "-")

    def normalize_button_identifier_keep_capitalization(self, button_name: str) -> str:
        """
        Takes BDD spec file input and converts it to the ID formatting in the data portal.
        Does not change the capitalization of the string.
        """
        return button_name.replace(" ", "-")

    def normalize_identifier_underscore(self, identifier_name: str) -> str:
        """Takes BDD spec file input and converts it to the ID formatting in the data portal"""
        return identifier_name.lower().replace(" ", "_")

    def normalize_identifier_underscore_keep_capitalization(self, identifier_name: str) -> str:
        """
        Takes BDD spec file input and converts it to the ID formatting in the data portal
        Does not change the capitalization of the string.
        """
        return identifier_name.replace(" ", "_")

    def normalize_applied_filter_name(self, filter_name: str) -> List[str]:
        periods = [char for char in filter_name if char == "."]
        filter_name = filter_name.replace("_", " ").replace(".", " ")
        filter_name = filter_name.split(" ")
        if len(periods) > 1:
            filter_name = filter_name[2:]
        filter_name = " ".join(word[0].upper() + word[1:] for word in filter_name)
        return filter_name

    def make_input_0_index(self, bdd_input):
        """Takes BDD spec file input and converts it to index 0"""
        bdd_input_int = int(bdd_input)
        bdd_input_int -= 1
        bdd_input_string = str(bdd_input_int)
        return bdd_input_string

    def strip_string_for_comparison(self, string_to_strip):
        """Takes a string and strips it down for comparison"""
        # Some create filtered cohort buttons have a denominator,
        # and we are only interested in the numerator for comparison.
        string_to_strip = string_to_strip.split("/")[0]
        string_to_strip = string_to_strip.replace(",", "")
        string_to_strip = string_to_strip.replace(" ", "")
        # Removes any non-numeral character
        string_to_strip = re.sub(r'\D', '', string_to_strip)
        return string_to_strip

    def get_text_by_data_testid(self, text_id_to_collect):
        """Returns text from given data-testid"""
        text_id_to_collect = self.normalize_button_identifier(text_id_to_collect)
        locator = GenericLocators.DATA_TESTID_TEXT_IDENT(text_id_to_collect)
        return self.get_text(locator)

    def get_cohort_bar_case_count(self):
        """Returns the count of cases in the current cohort"""
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        cohort_bar_case_count = self.get_text(
            GenericLocators.TEXT_COHORT_BAR_CASE_COUNT
        )
        return cohort_bar_case_count

    def get_showing_count_text(self):
        """Returns the text of how many items are being shown on the page"""
        locator = GenericLocators.SHOWING_NUMBER_OF_ITEMS
        return self.get_text(locator)

    def get_table_showing_count_text(self, table_name):
        """Returns the text of how many items are being shown on the specified table"""
        table_name = self.normalize_button_identifier(table_name)
        locator = GenericLocators.SHOWING_NUMBER_OF_ITEMS_IN_TABLE(table_name)
        return self.get_text(locator)

    def get_table_item_count_text(self, table_name):
        """Returns how many items are in the table by retrieving the number in the 'Total Of xx Object' string"""
        table_name = self.normalize_button_identifier(table_name)
        locator = GenericLocators.TABLE_ITEM_COUNT_IDENT(table_name)

        # When this number loads, it will display 0 before showing the actual count.
        # So, we check to see if the text displays '0' or '--' which means it's still loading.
        # We check for the text to display an actual number for approximately 5 seconds.
        # If no number appears, we return the '0' as is. In rare scenarios, the count
        # really is 0, so this does not necessarily mean failure.
        item_count = self.get_text(locator)
        retry_counter = 0
        while (item_count == "0" or item_count =="--"):
            time.sleep(1)
            item_count = self.get_text(locator)
            retry_counter = retry_counter+1
            if retry_counter >= 5:
                break
        return self.get_text(locator)

    def get_filter_selection_count(self, filter_group_name, selection):
        """Returns the count of how many items are associated with that filter in the current cohort"""
        locator = GenericLocators.FILTER_GROUP_SELECTION_COUNT_IDENT(
            filter_group_name, selection
        )
        return self.get_text(locator)

    def get_table_header_text_by_column(self, column):
        """
        Gets text of table header by column.
        Column indexing begins at '1'
        """
        table_header_text_locator = GenericLocators.TEXT_TABLE_HEADER(column)
        return self.get_text(table_header_text_locator)

    def get_table_header_text_by_column_in_specified_table(self, table_name, column):
        """
        Gets text of specified table header by column.
        Column indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_header_text_locator = GenericLocators.TEXT_SPECIFIED_TABLE_HEADER(table_name,column)
        return self.get_text(table_header_text_locator)

    def get_table_header_text_by_row_in_specified_vertical_table(self, table_name, column):
        """
        Gets text of specified table header by row. For tables that are vertically aligned.
        Row indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_header_text_locator = GenericLocators.TEXT_SPECIFIED_VERTICAL_TABLE_HEADER(table_name,column)
        return self.get_text(table_header_text_locator)

    def get_table_body_text_by_row_column(self, row, column):
        """
        Gets text from table body by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_SELECT(row, column)
        return self.get_text(table_locator_to_select)

    def get_specified_table_body_text_by_row_column(self, table_name, row, column):
        """
        In specified table, gets text from table body by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_SELECT_IN_SPECIFIED_TABLE(table_name, row, column)
        # Try to return drilled down locator.
        if self.is_visible(table_locator_to_select):
            return self.get_inner_text(table_locator_to_select)
        # If that is not available, return a higher level locator.
        else:
            table_locator_to_select = GenericLocators.TABLE_AREA_TO_COLLECT_IN_SPECIFIED_TABLE(table_name, row, column)
            return self.get_inner_text(table_locator_to_select)

    def hover_table_body_by_row_column(self, row, column):
        """
        Hovers over specified cell in table by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_SELECT(row, column)
        # The hover function is finicky. Calling it twice, and then bypassing
        # actionability checks seem to make the hover more consistent.

        self.hover(table_locator_to_select)
        self.hover(table_locator_to_select, force=True)

    def hover_table_body_by_row_column_in_specified_table(self, table_name, row, column):
        """
        In specified table, hovers over specified cell in table by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_SELECT_IN_SPECIFIED_TABLE(table_name, row, column)
        # The hover function is finicky. Calling it twice, and then bypassing
        # actionability checks seem to make the hover more consistent.

        self.hover(table_locator_to_select)
        self.hover(table_locator_to_select, force=True)

    def wait_until_locator_is_visible(self, locator, timeout=30000):
        """wait for element to have non-empty bounding box and no visibility:hidden"""
        self.driver.locator(locator).wait_for(state="visible", timeout=timeout)

    def wait_until_locator_is_detached(self, locator):
        """wait for element to not be present in DOM"""
        self.driver.locator(locator).wait_for(state="detached", timeout=60000)

    def wait_until_locator_is_hidden(self, locator):
        """wait for element to be either detached from DOM, or have an empty bounding box or visibility:hidden"""
        self.driver.locator(locator).wait_for(state="hidden", timeout=15000)

    def wait_for_text_in_temporary_message(self, text, action="remove modal"):
        """
        Waits for a piece of text to appear in the temporary cohort modal.
        That modal appears after an action has been performed on a cohort
        state (e.g create, save, delete, etc).
        """
        text_locator = GenericLocators.TEXT_IN_PARAGRAPH(text)
        try:
            self.wait_until_locator_is_visible(text_locator)
            if action.lower() == "remove modal":
                # On occasion, the automation will move so fast and click the close 'x' button
                # it changes what the active cohort is. I cannot reproduce it manually, and it stops
                # when I put this sleep here.
                time.sleep(1)
                x_button_locator = GenericLocators.BUTTON_CLOSE_NOTIFICATION
                # Remove the message after locating it.
                # The messages can pile up, so removing them is sometimes necessary for subsequent scenarios
                self.click(x_button_locator)
        except:
            return False
        return True

    def wait_for_selector(self, locator, timeout=30000):
        self.driver.wait_for_selector(locator, timeout=timeout)

    def wait_for_loading_spinner_to_be_visible(self, timeout=30000):
        locator = GenericLocators.LOADING_SPINNER_GENERIC
        self.wait_until_locator_is_visible(locator, timeout)

    def wait_for_loading_spinner_to_detatch(self):
        """Waits for the generic loading spinner to disappear on the page"""
        locator = GenericLocators.LOADING_SPINNER_GENERIC
        self.wait_until_locator_is_detached(locator)

    def wait_for_loading_spinner_cohort_bar_case_count_to_detatch(self):
        """Waits for the cohort bar case count loading spinner to disappear on the page"""
        locator = GenericLocators.LOADING_SPINNER_COHORT_BAR_CASE_COUNT
        self.wait_until_locator_is_detached(locator)

    def wait_for_loading_spinner_table_to_detatch(self):
        """Waits for the table (repository, projects, mutation frequency, etc.) loading spinner to disappear on the page"""
        locator = GenericLocators.LOADING_SPINNER_TABLE
        self.wait_until_locator_is_detached(locator)

    def wait_for_loading_spinners_to_detach(self):
        """
        We often have to wait for many possible loading spinners to detach.
        This function is a convenient way to do that at once.
        """
        time.sleep(0.5)
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_table_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_table_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()

    def wait_for_data_testid_to_be_visible(self, locator):
        """Normalizes a data-testid and waits for it to be visible"""
        normalized_locator = self.normalize_button_identifier(locator)
        locator = GenericLocators.DATA_TEST_ID_IDENT(normalized_locator)
        try:
            self.wait_until_locator_is_visible(locator, 45000)
        except:
            return False
        return True

    def wait_for_table_body_text_by_row_column(self, text, row, column):
        """
        Waits for text from table body by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_locator_to_select = GenericLocators.TABLE_TEXT_TO_WAIT_FOR(
            text, row, column
        )
        self.wait_until_locator_is_visible(table_locator_to_select)

    def wait_for_specified_table_body_text_by_row_column(self, table_name, text, row, column):
        """
        In specified table, waits for text from table body by giving a row and column.
        Row and Column indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_locator_to_select = GenericLocators.TABLE_TEXT_TO_WAIT_FOR_IN_SPECIFIED_TABLE(
            table_name, text, row, column
        )
        self.wait_until_locator_is_visible(table_locator_to_select)

    def is_text_present(self, text):
        locator = GenericLocators.TEXT_IDENT(text)
        try:
            self.wait_until_locator_is_visible(locator)
        except:
            return False
        return True

    def is_text_not_present(self, text):
        locator = GenericLocators.TEXT_IDENT(text)
        try:
            self.wait_until_locator_is_hidden(locator)
        except:
            return False
        return True

    def is_message_id_text_present(self, message_id, text):
        message_id = self.normalize_button_identifier(message_id)
        locator = GenericLocators.DATA_TESTID_TEXT_IDENT(message_id)
        try:
            self.wait_until_locator_is_visible(locator)
            message_text = self.get_text(locator)
            message_text == text
        except:
            return False
        return True

    def is_cohort_bar_case_count_present(self, case_count):
        locator = GenericLocators.COHORT_BAR_CASE_COUNT(case_count)
        try:
            self.wait_until_locator_is_visible(locator)
        except:
            return False
        return True

    def is_data_testid_present(self, data_testid):
        locator = GenericLocators.DATA_TEST_ID_IDENT(data_testid)
        is_data_testid_present = self.is_visible(locator)
        return is_data_testid_present

    def is_data_testid_button_present(self, data_testid):
        locator = GenericLocators.DATA_TESTID_BUTTON_IDENT(data_testid)
        is_data_testid_present = self.is_visible(locator)
        return is_data_testid_present

    def is_button_disabled(self, button_name):
        """Returns if the data-testid button is disabled"""
        button_name = self.normalize_button_identifier(button_name)
        locator = GenericLocators.DATA_TESTID_BUTTON_IDENT(button_name)
        is_button_disabled = self.is_disabled(locator)
        return is_button_disabled

    def is_named_button_disabled(self, button_name):
        """Returns if the specified named button is disabled"""
        locator = GenericLocators.BUTTON_BY_DISPLAYED_TEXT(button_name)
        is_button_disabled = self.is_disabled(locator)
        return is_button_disabled

    def is_dropdown_option_text_button_disabled(self, dropdown_text_option):
        locator = GenericLocators.BUTTON_TEXT_DROPDOWN_MENU_OPTION(dropdown_text_option)
        is_button_disabled = self.is_disabled(locator)
        return is_button_disabled

    def is_button_area_expanded(self, button_name):
        """Returns if the data-testid button is expanded"""
        button_name = self.normalize_button_identifier(button_name)
        locator = GenericLocators.DATA_TESTID_BUTTON_IDENT(button_name)
        return self.get_attribute(locator,"aria-expanded")

    def is_cart_count_correct(self, correct_file_count):
        """Returns if cart count is correct"""
        locator = GenericLocators.CART_IDENT
        cart_text = self.get_text(locator)
        cart_count = cart_text.replace("Cart", "")
        return cart_count == correct_file_count

    # Checks to see if specified filter card is present
    def is_filter_card_present(self, filter_group_name):
        locator = GenericLocators.FILTER_GROUP_IDENT(filter_group_name)
        try:
            self.wait_until_locator_is_visible(locator)
        except:
            return False
        return True

    def is_show_more_or_show_less_button_visible_within_filter_card(self, facet_group_name, label):
        """Returns if the 'show more or show less' button is visible on a facet card"""
        locator = GenericLocators.FILTER_GROUP_SHOW_MORE_LESS_IDENT(
            facet_group_name, label
        )
        return self.is_visible(locator)

    def is_custom_filter_in_properties_table_present(self, button_name):
        """Returns if custom filter is present in properties table"""
        locator = GenericLocators.BUTTON_CUSTOM_FILTER_IN_PROPERTIES_TABLE(button_name)
        return self.is_visible(locator)

    def is_filter_card_custom_filter_text_present(self, facet_card, text):
        """Returns if filter text is present on given facet card"""
        locator = GenericLocators.FILTER_GROUP_CUSTOM_FILTER_TEXT_IDENT(facet_card, text)
        result = self.is_visible(locator)
        return result

    def is_no_active_cohort_filter_text_present(self):
        """
        Returns if the text 'No filters currently applied.' is displayed
        in the active cohort filter area
        """
        text_no_active_filter_locator = GenericLocators.TEXT_NO_ACTIVE_COHORT_FILTERS
        return self.is_visible(text_no_active_filter_locator)

    def is_table_displaying_text(self, table_id, table_text):
        table_id = self.normalize_button_identifier(table_id)
        table_text_locator = GenericLocators.TABLE_TEXT_IDENT(table_id, table_text)
        try:
            self.wait_until_locator_is_visible(table_text_locator)
        except:
            return False
        return True

    def is_table_displaying_text_no_wait(self, table_id, table_text):
        """Returns if the given table is displaying specified text. This does not wait
        for text to appear, it returns immediately if is displayed or not."""
        table_id = self.normalize_button_identifier(table_id)
        table_text_locator = GenericLocators.TABLE_TEXT_IDENT(table_id, table_text)
        return self.is_visible(table_text_locator)

    def is_temporary_message_present(self):
        """Returns if temporary message modal is present"""
        locator = GenericLocators.BUTTON_CLOSE_NOTIFICATION
        return self.is_visible(locator)

    def click_data_testid(self, data_testid):
        locator = GenericLocators.DATA_TEST_ID_IDENT(data_testid)
        self.click(locator)

    def click_button_data_testid(self, data_testid):
        locator = GenericLocators.DATA_TESTID_BUTTON_IDENT(data_testid)
        self.click(locator)

    def click_button_data_testid_normalize(self, data_testid):
        data_testid = self.normalize_button_identifier(data_testid)
        locator = GenericLocators.DATA_TESTID_BUTTON_IDENT(data_testid)
        self.click(locator)

    def click_button_with_displayed_text_name(self, button_text_name):
        """Clicks a button based on its displayed text"""
        locator = GenericLocators.BUTTON_BY_DISPLAYED_TEXT(button_text_name)
        self.click(locator)

    def click_button_in_modal_with_displayed_text_name(self, button_text_name):
        """Clicks a button in a modal based on its displayed text"""
        locator = GenericLocators.BUTTON_IN_MODAL_BY_DISPLAYED_TEXT(button_text_name)
        self.click(locator)

    def click_button_ident_a_with_displayed_text_name(self, button_text_name):
        """Clicks a button based on its displayed text with a CSS tag of 'a'"""
        locator = GenericLocators.BUTTON_A_BY_TEXT_IDENT(button_text_name)
        self.click(locator)

    def click_footer_button_ident_with_displayed_text_name(self, button_text_name):
        """Clicks a button in the footer based on its displayed text"""
        locator = GenericLocators.BUTTON_IN_FOOTER_BY_TEXT_IDENT(button_text_name)
        self.click(locator)

    def click_link_data_testid(self, link_data_testid):
        """Clicks a link with a data-testid"""
        link_data_testid = self.normalize_button_identifier(link_data_testid)
        locator = GenericLocators.DATA_TESTID_LINK_IDENT(link_data_testid)
        self.click(locator)

    def click_checkbox(self, checkbox_name):
        """Clicks a checkbox with given, unmodified name"""
        locator = GenericLocators.CHECKBOX_IDENT(checkbox_name)
        self.click(locator)

    def click_checkbox_normalize_ident(self, checkbox_name):
        """Clicks a checkbox with given, unmodified name"""
        checkbox_name = self.normalize_button_identifier(checkbox_name)
        locator = GenericLocators.CHECKBOX_IDENT(checkbox_name)
        self.click(locator)

    def click_radio_button(self, radio_name):
        """Clicks a radio button in a filter card"""
        locator = GenericLocators.RADIO_BUTTON_IDENT(radio_name)
        self.click(locator)

    def click_radio_data_testid_button(self, radio_name):
        """Clicks a data-testid radio button"""
        locator = GenericLocators.RADIO_DATA_TESTID_BUTTON_IDENT(radio_name)
        self.click(locator)

    def click_close_modal_button(self):
        """Clicks 'X' to close a modal"""
        locator = GenericLocators.BUTTON_CLOSE_MODAL
        self.click(locator)

    def click_close_temporary_message(self):
        """Clicks 'X' to close a modal"""
        locator = GenericLocators.BUTTON_CLOSE_NOTIFICATION
        self.click(locator)
        time.sleep(0.5)

    def click_undo_in_message(self):
        """Clicks 'undo' in a modal message"""
        locator = GenericLocators.UNDO_BUTTON_IN_TEMP_MESSAGE
        self.click(locator)

    def click_set_as_current_cohort_in_message(self):
        """Clicks 'Set this as your current cohort' in a modal message"""
        locator = GenericLocators.SET_AS_CURRENT_COHORT_IN_TEMP_MESSAGE
        self.click(locator)

    def click_create_or_save_button_in_cohort_modal(self):
        """Clicks 'Create' or 'Save' in cohort modal"""
        locator = GenericLocators.BUTTON_CREATE_OR_SAVE_COHORT_IN_MODAL
        self.click(locator)

    def click_text_option_from_dropdown_menu(self, dropdown_option):
        """Clicks a text option from a dropdown menu"""
        locator = GenericLocators.TEXT_DROPDOWN_MENU_OPTION(dropdown_option)
        self.click(locator)

    def click_has_text_option_from_dropdown_menu(self, dropdown_option):
        """Clicks a text option from a dropdown menu based off given partial text"""
        dropdown_approximate_locator = GenericLocators.TEXT_HAS_DROPDOWN_MENU_OPTION(dropdown_option)
        self.click(dropdown_approximate_locator)

    def click_button_dropdown_in_graph(self, graph_name):
        """Clicks dropdown button in specified graph"""
        graph_name = self.normalize_button_identifier(graph_name)
        locator = GenericLocators.BUTTON_IN_GRAPH(graph_name)
        self.click(locator)

    def click_button_in_most_frequent_somatic_mutations(self, button_name):
        """Clicks download button in most frequent somatic mutations table"""
        locator = GenericLocators.BUTTON_IN_MOST_FREQUENT_SOMATIC_MUTATIONS_TABLE(button_name)
        self.click(locator)

    def click_text_option_from_dropdown_data_image(self, dropdown_option):
        """Clicks a text option from a graph dropdown menu"""
        locator = GenericLocators.BUTTON_TEXT_DROPDOWN_IMAGE_DATA(dropdown_option)
        self.click(locator)

    def clear_active_cohort_filters(self):
        """
        Clears the active cohort filters by clicking the "Clear All" button
        """
        if not self.is_no_active_cohort_filter_text_present():
            button_clear_all_active_cohort_filters_locator = (
                GenericLocators.BUTTON_CLEAR_ACTIVE_COHORT_FILTERS
            )
            self.click(button_clear_all_active_cohort_filters_locator)
            text_no_active_filter_locator = (
                GenericLocators.TEXT_NO_ACTIVE_COHORT_FILTERS
            )
            self.wait_for_data_testid_to_be_visible(text_no_active_filter_locator)
            self.wait_for_loading_spinners_to_detach()

    def click_column_selector_button(self):
        """Clicks table column selector button"""
        locator = GenericLocators.BUTTON_COLUMN_SELECTOR
        self.click(locator)

    def click_column_selector_button_in_specified_table(self, table_name):
        """In specified table, clicks table column selector button"""
        table_name = self.normalize_button_identifier(table_name)
        locator = GenericLocators.BUTTON_COLUMN_SELECTOR_IN_TABLE(table_name)
        self.click(locator)

    def click_switch_for_column_selector(self, switch_name):
        """In the column selector pop-up modal, clicks specified switch"""
        switch_name = self.normalize_identifier_underscore(switch_name)
        locator = GenericLocators.SWITCH_COLUMN_SELECTOR(switch_name)
        self.click(locator)

    def click_reset_column_select_options(self):
        locator = GenericLocators.BUTTON_RESET_COLUMN_SELECTIONS
        self.click(locator)

    def change_number_of_entries_shown(self, entries_to_show):
        """
        Changes number of entries shown in the table using the show entries button,
        and selecting an option from the dropdown list.
        """
        entries_button_locator = GenericLocators.BUTTON_ENTRIES_SHOWN
        self.click(entries_button_locator)

        dropdown_entries_to_show_locator = (
            GenericLocators.DROPDOWN_LIST_CHANGE_NUMBER_OF_ENTRIES_SHOWN(
                entries_to_show
            )
        )
        self.click(dropdown_entries_to_show_locator)

    def change_number_of_entries_shown_in_specified_table(self, table_name, entries_to_show):
        """
        Changes number of entries shown in the table using the show entries button,
        and selecting an option from the dropdown list.
        """
        table_name = self.normalize_button_identifier(table_name)
        entries_button_locator = GenericLocators.BUTTON_ENTRIES_SHOWN_IN_SPECIFIED_TABLE(table_name)
        self.click(entries_button_locator)

        dropdown_entries_to_show_locator = (
            GenericLocators.DROPDOWN_LIST_CHANGE_NUMBER_OF_ENTRIES_SHOWN_IN_SPECIFIED_TABLE(
                table_name, entries_to_show
            )
        )
        self.click(dropdown_entries_to_show_locator)

    def make_selection_within_filter_group(self, filter_group_name, selection):
        """Clicks a checkbox within a filter group"""
        locator = GenericLocators.FILTER_GROUP_SELECTION_IDENT(
            filter_group_name, selection
        )
        self.click(locator)

    def perform_action_within_filter_card(self, filter_group_name, action):
        """Performs an action in a filter group e.g sorting, resetting, flipping the chart, etc."""
        locator = GenericLocators.FILTER_GROUP_ACTION_IDENT(filter_group_name, action)
        self.click(locator)

    def flip_switch_in_filter_card(self, filter_group_name):
        """Flips the switch in specified filter card"""
        locator = GenericLocators.FILTER_GROUP_FLIP_SWITCH_IDENT(filter_group_name)
        self.click(locator)

    def click_show_more_less_within_filter_card(self, filter_group_name, label):
        """Clicks the show more or show less object"""
        locator = GenericLocators.FILTER_GROUP_SHOW_MORE_LESS_IDENT(
            filter_group_name, label
        )
        self.click(locator)

    def type_in_filter_card_search_text_area(self, facet_group_name, textbox_id, text):
        """Send keys in the search textbox area"""
        textbox_id = self.normalize_button_identifier(textbox_id)
        locator = GenericLocators.FILTER_GROUP_TEXT_AREA_IDENT(
            facet_group_name, textbox_id
        )
        self.send_keys(locator, text)

    def select_table_by_row_column(self, row, column):
        """
        Selects values from tables by giving a row and column
        Row and Column indexing begins at '1'
        """
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_CLICK(row, column)
        self.hover(table_locator_to_select)
        self.click(table_locator_to_select, True)

    def select_specified_table_by_row_column(self, table_name, row, column):
        """
        In specified table, selects values from tables by giving a row and column
        Row and Column indexing begins at '1'
        """
        table_name = self.normalize_button_identifier(table_name)
        table_checkbox_to_click = GenericLocators.TABLE_CHECKBOX_TO_CLICK_IN_SPECIFIED_TABLE(table_name, row, column)
        table_locator_to_select = GenericLocators.TABLE_AREA_TO_CLICK_IN_SPECIFIED_TABLE(table_name, row, column)
        # Try to click potential checkbox first.
        if self.is_visible(table_checkbox_to_click):
            self.hover(table_checkbox_to_click)
            self.click(table_checkbox_to_click, True)
        # Try to click drilled down locator.
        elif self.is_visible(table_locator_to_select):
            self.hover(table_locator_to_select)
            self.click(table_locator_to_select, True)
        # If that is not available, click a higher level locator.
        else:
            table_locator_to_select = GenericLocators.TABLE_AREA_TO_SELECT_IN_SPECIFIED_TABLE(table_name, row, column)
            self.hover(table_locator_to_select)
            self.click(table_locator_to_select, True)

    def send_text_into_search_bar(self, text_to_send, aria_label):
        """Sends text into search bar based on its aria_label"""
        locator = GenericLocators.SEARCH_BAR_ARIA_IDENT(aria_label)
        self.wait_until_locator_is_visible(locator)
        self.send_keys(locator, text_to_send)

    def send_text_into_text_box(self, text_to_send, text_box_id):
        """Sends text into data-testid textbox"""
        text_box_id = self.normalize_button_identifier(text_box_id)
        locator = GenericLocators.TEXT_BOX_IDENT(text_box_id)
        self.wait_until_locator_is_visible(locator)
        self.send_keys(locator, text_to_send)

    def send_text_into_specified_table_search_bar(self, table_specified, text_to_send):
        """Sends text into a specified table search bar"""
        table_specified = self.normalize_button_identifier(table_specified)
        locator = GenericLocators.SEARCH_BAR_IN_SPECIFIED_TABLE_IDENT(table_specified)
        self.wait_until_locator_is_visible(locator)
        self.send_keys(locator, text_to_send)

    def send_text_into_table_search_bar(self, text_to_send):
        """Sends text into a table search bar"""
        locator = GenericLocators.SEARCH_BAR_TABLE_IDENT
        self.wait_until_locator_is_visible(locator)
        self.send_keys(locator, text_to_send)

    def quick_search_and_click(self, text):
        """
        Sends text into the quick search bar in the upper-right corner of the data portal.
        Then clicks the first result in the search result area. Best used with a UUID.
        """
        self.send_keys(GenericLocators.QUICK_SEARCH_BAR_IDENT, text)
        time.sleep(1)
        first_result_locator = GenericLocators.QUICK_SEARCH_BAR_FIRST_RESULT
        self.click(first_result_locator)

    def global_quick_search(self, text):
        """Sends text into the quick search bar in the upper-right corner of the data portal."""
        self.send_keys(GenericLocators.QUICK_SEARCH_BAR_IDENT, text)

    def validate_global_quick_search_result_text(self, result_in_list, text):
        """Specifies a result from the quick search bar result list. Validates expected text is present."""
        result_in_list = self.make_input_0_index(result_in_list)
        locator_result_category = GenericLocators.QUICK_SEARCH_BAR_RESULT_TEXT(
            result_in_list, text
        )
        self.wait_until_locator_is_visible(locator_result_category)

    def click_global_quick_search_result(self, result_in_list):
        """Specifies a result from the quick search bar result list. Clicks that result."""
        result_in_list = self.make_input_0_index(result_in_list)
        locator_result_in_list = GenericLocators.QUICK_SEARCH_BAR_NUMBERED_RESULT(
            result_in_list
        )
        self.click(locator_result_in_list)

    # This section of functions is for handling new tabs
    def perform_action_handle_new_tab(self, source: str, button: str):
        """
        perform_action_handle_new_tab performs an action to open a new tab,
        and then returns a page object for that new tab.

        :param source: Specifies what function is causing the action to open the new tab
        :param button: ID or Name of the button that is being clicked to open the new tab
        :return: a page object for the new tab that has been opened
        """
        sources = {
            "Home Page": self.click_button_ident_a_with_displayed_text_name,
            "Footer": self.click_footer_button_ident_with_displayed_text_name,
            "Cart": self.click_link_data_testid,
        }
        driver = WebDriver.page
        with driver.context.expect_page() as tab:
            sources.get(source)(button)
        new_tab = tab.value
        return new_tab

    def click_in_table_handle_new_tab(self, table_name, row, column):
        """
        click_in_table_handle_new_tab clicks in a table to open a new tab,
        and then returns a page object for that new tab.

        :param table_name: Specifies what table to click on
        :param row: Row number to click
        :param column: Column number to click
        :return: a page object for the new tab that has been opened
        """
        driver = WebDriver.page
        with driver.context.expect_page() as tab:
            self.select_specified_table_by_row_column(table_name, row, column)
        new_tab = tab.value
        return new_tab

    def is_text_visible_on_new_tab(self, new_tab, text_to_check):
        """
        is_text_visible_on_new_tab checks for text on a given tab page.

        :param new_tab: The tab page to be checked.
        :param text_to_check: The <p> text to be searched for.
        """
        expected_text_locator = GenericLocators.TEXT_IN_PARAGRAPH(text_to_check)
        try:
            new_tab.locator(expected_text_locator).wait_for(state="visible")
        except:
            return False
        return True
