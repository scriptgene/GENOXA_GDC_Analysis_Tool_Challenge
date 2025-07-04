from playwright.sync_api import Page

from ....base.base_page import BasePage

import time

class HeaderSectionLocators:
    BUTTON_IDENT = lambda button_name: f"[data-testid='button-header-{button_name}']"
    BUTTON_LOGIN = "[data-testid='button-header-login'] >> nth=1"
    BUTTON_USERNAME = "[data-testid='button-header-username'] >> nth=1"

    # These are all locators that only appear when the respective page has fully loaded
    ANALYSIS_CENTER_WAIT_FOR_ELEMENT = "[data-testid='Clinical Data Analysis-tool']"
    PROJECTS_WAIT_FOR_ELEMENT = "[data-testid='button-json-projects-table']"
    COHORT_BUILDER_WAIT_FOR_ELEMENT = "[data-testid='button-cohort-builder-general']"
    COHORT_BUILDER_ADDITIONAL_WAIT_FOR_ELEMENT = "[data-testid='title-cohort-builder-facet-groups']"
    REPOSITORY_WAIT_FOR_ELEMENT = "[data-testid='button-json-files-table']"
    REPOSITORY_ADDITIONAL_WAIT_FOR_ELEMENT = "[data-testid='text-showing-count']"
    HOME_WAIT_FOR_ELEMENT = "[data-testid='text-cases-gdc-count']"
    MANAGE_SETS_WAIT_FOR_ELEMENT = "[data-testid='button-create-set']"
    BROWSE_ANNOTATIONS_WAIT_FOR_ELEMENT = "[data-testid='table-annotations']"


class HeaderSection(BasePage):
    def __init__(self, driver: Page, url):
        self.driver = driver  # driver is PW page

    def navigate_to_main_pages(self, button_name: str):
        button_name = self.normalize_button_identifier(button_name)
        locator = HeaderSectionLocators.BUTTON_IDENT(button_name)
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        self.wait_until_locator_is_visible(locator)
        self.click(locator, True, 15000)
        try:
            self.wait_for_page_to_load(button_name)
        except:
            self.click(locator, True, 15000)
            self.wait_for_page_to_load(button_name)

    # Pages in the data portal do not load instantaneously.
    # We want to wait for the main content of the page to load before continuing the test.
    def wait_for_page_to_load(self, page_to_load):
        page_to_load = page_to_load.lower()
        if page_to_load == "analysis":
            self.wait_for_selector(
                HeaderSectionLocators.ANALYSIS_CENTER_WAIT_FOR_ELEMENT
            )
            # Wait for the central, initial loading spinner that appears.
            # If it doesn't appear, we move on.
            try:
                self.wait_for_loading_spinner_to_be_visible(5000)
                self.wait_for_loading_spinner_to_detatch()
            except:
                pass
        elif page_to_load == "projects":
            self.wait_for_selector(HeaderSectionLocators.PROJECTS_WAIT_FOR_ELEMENT)
        elif page_to_load == "cohort":
            self.wait_for_selector(
                HeaderSectionLocators.COHORT_BUILDER_WAIT_FOR_ELEMENT
            )
            self.wait_for_selector(
                HeaderSectionLocators.COHORT_BUILDER_ADDITIONAL_WAIT_FOR_ELEMENT
            )
        elif page_to_load == "downloads":
            # Repository page does not load quickly, and automation will move too fast at times
            self.wait_for_selector(HeaderSectionLocators.REPOSITORY_WAIT_FOR_ELEMENT)
            self.wait_for_selector(
                HeaderSectionLocators.REPOSITORY_ADDITIONAL_WAIT_FOR_ELEMENT
            )
        elif page_to_load == "home":
            self.wait_for_selector(HeaderSectionLocators.HOME_WAIT_FOR_ELEMENT)
        elif page_to_load == "manage-sets":
            self.wait_for_selector(HeaderSectionLocators.MANAGE_SETS_WAIT_FOR_ELEMENT)
        elif page_to_load == "cart":
            # We are not guaranteed to see any one thing on the cart page, so instead
            # we will pause for a moment.
            time.sleep(1)
        elif page_to_load == "browse-annotations":
            self.wait_for_selector(HeaderSectionLocators.BROWSE_ANNOTATIONS_WAIT_FOR_ELEMENT)

        self.wait_for_loading_spinner_table_to_detatch()
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        self.wait_for_loading_spinner_table_to_detatch()
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()

    def login_to_data_portal_if_possible(self):
        login_locator = HeaderSectionLocators.BUTTON_LOGIN
        if self.is_visible(login_locator):
            self.click(login_locator)
            wait_for_username_locator = HeaderSectionLocators.BUTTON_USERNAME
            self.wait_until_locator_is_visible(wait_for_username_locator, 60000)
            self.wait_for_loading_spinners_to_detach()
