import time

from playwright.sync_api import Page

from ....base.base_page import BasePage
from step_impl.apps.gdc_data_portal_v2.pages.header_section import HeaderSectionLocators
from step_impl.apps.gdc_data_portal_v2.pages.cohort_comparison_page import (
    CohortComparisonLocators,
)


class AnalysisCenterLocators:
    BUTTON_APP_PLAY_OR_DEMO = lambda app_name: f'[data-testid="button-{app_name}"]'

    SELECT_DESCRIPTION_TOOL = (
        lambda tool_name: f'[data-testid="{tool_name}-tool"] >> [data-testid="select-description-tool"]'
    )
    TEXT_DESCRIPTION_TOOL = (
        lambda tool_name: f'[data-testid="{tool_name}-tool"] >> [data-testid="text-description-tool"]'
    )
    TEXT_CASES_COUNT_ON_TOOL_CARD = (
        lambda tool_name: f'[data-testid="{tool_name}-tool"] >> [data-testid="text-case-count-tool"]'
    )
    TOOLTIP_ZERO_CASES_ON_TOOL_CARD = (
        lambda tool_name: f'[data-testid="{tool_name}-tool"] [data-testid="text-case-count-tool"] svg'
    )

    CORE_TOOL_PROJECTS = '[data-testid="button-core-tools-Projects"]'
    CORE_TOOL_COHORT_BUILDER = '[data-testid="button-core-tools-Cohort Builder"]'
    CORE_TOOL_REPOSITORY = '[data-testid="button-core-tools-Repository"]'

    ANALYSIS_CENTER_HEADER = '[data-testid="button-header-analysis"]'

    MUTATION_FREQUENCY_WAIT_FOR_ELEMENT = "[data-testid='button-mutations-tab']"
    CLINICAL_DATA_ANALYSIS_WAIT_FOR_ELEMENT = "[data-testid='Gender-card']"


class AnalysisCenterPage(BasePage):
    def __init__(self, driver: Page, url: str) -> None:
        self.URL = "{}/analysis_page".format(url)
        self.driver = driver  # driver is PW page

    def visit(self):
        self.driver.goto(self.URL)

    def is_analysis_center_page_present(self):
        locator = AnalysisCenterLocators.CORE_TOOL_REPOSITORY
        return self.is_visible(locator)

    def navigate_to_app(self, app_name: str):
        locator = AnalysisCenterLocators.BUTTON_APP_PLAY_OR_DEMO(app_name)
        self.click(locator)
        self.wait_for_app_page_to_load(app_name)

    def click_analysis_tool_description(self, tool_name):
        """Clicks the description for an analysis tool"""
        locator = AnalysisCenterLocators.SELECT_DESCRIPTION_TOOL(tool_name)
        self.click(locator)

    def get_analysis_tool_text(self, tool_name):
        """After an analysis tool description is made visible, this returns the description"""
        locator = AnalysisCenterLocators.TEXT_DESCRIPTION_TOOL(tool_name)
        return self.get_text(locator)

    def get_analysis_tool_cases_count(self, tool_name):
        """Returns case count on given analysis tool"""
        self.wait_for_loading_spinner_to_detatch()
        locator = AnalysisCenterLocators.TEXT_CASES_COUNT_ON_TOOL_CARD(tool_name)
        # If we get " " or " Cases" on return, that means the analysis card is still loading information. We wait until
        # it fully loads or 10 seconds elapses before returning text.
        retry_counter = 0
        while ((self.get_text(locator) == " ") or (self.get_text(locator) == " Cases")):
            time.sleep(1)
            retry_counter = retry_counter+1
            if retry_counter >= 10:
                break
        cases_count = self.get_text(locator)
        # Remove the "Cases" part of the string. We do not need that for comparison.
        cases_count = cases_count.replace("Cases", "")
        cases_count = cases_count.replace(" ", "")
        return cases_count

    def get_analysis_tool_tooltip_text(self, tool_name, tooltip_text):
        """
        When an analysis tool has 0 cases, a tooltip appears on the tool card.
        This hovers over the tooltip and checks if the text matches our spec file
        """
        locator = AnalysisCenterLocators.TOOLTIP_ZERO_CASES_ON_TOOL_CARD(tool_name)
        # The click and sleep is to make the previous tooltip disappear, so it does not
        # affect the next tooltip check
        self.click(locator)
        time.sleep(0.5)
        self.hover(locator)
        is_tooltip_text_present = self.is_text_present(tooltip_text)
        return is_tooltip_text_present

    def core_tools_navigation_check(self):
        # First element in the set: a featured tool navigate button
        # Second element in the set: an element to check if user landed on correct page
        nav_and_location = [
            (
                AnalysisCenterLocators.CORE_TOOL_PROJECTS,
                HeaderSectionLocators.PROJECTS_WAIT_FOR_ELEMENT,
            ),
            (
                AnalysisCenterLocators.CORE_TOOL_COHORT_BUILDER,
                HeaderSectionLocators.COHORT_BUILDER_WAIT_FOR_ELEMENT,
            ),
            (
                AnalysisCenterLocators.CORE_TOOL_REPOSITORY,
                HeaderSectionLocators.REPOSITORY_WAIT_FOR_ELEMENT,
            ),
        ]
        # Click on a featured tool nav button, then validate user arrived on correct page
        for navigation, location in nav_and_location:
            try:
                self.click(navigation)
                self.wait_until_locator_is_visible(location)
                # Navigate back to the analysis center for the next test
                self.click(AnalysisCenterLocators.ANALYSIS_CENTER_HEADER, True)
                self.wait_until_locator_is_visible(
                    HeaderSectionLocators.ANALYSIS_CENTER_WAIT_FOR_ELEMENT
                )
            except:
                return False
        return True

    # Pages in the data portal do not load instantaneously.
    # We want to wait for the main content of the page to load before continuing the test.
    def wait_for_app_page_to_load(self, page_to_load):
        page_to_load = page_to_load.lower()
        if (
            page_to_load == "mutation frequency"
            or page_to_load == "mutation frequency demo"
        ):
            # Been running into loading failures on mutation frequency.
            # Usually happens when running regression. The page will not load at all.
            # Simple solution I can think of is refreshing the page if that happens.
            try:
                self.wait_for_selector(
                    AnalysisCenterLocators.MUTATION_FREQUENCY_WAIT_FOR_ELEMENT, 15000
                )
            except:
                self.driver.reload()
                time.sleep(10)
                self.wait_for_selector(
                    AnalysisCenterLocators.MUTATION_FREQUENCY_WAIT_FOR_ELEMENT, 30000
                )

        if page_to_load == "cohort comparison demo":
            # Need to wait for loading spinners to be present, for them to disappear,
            # and wait for a special loading spinner attached to the survival plot to disappear
            try:
                self.wait_for_loading_spinner_to_be_visible(15000)
                self.wait_for_loading_spinner_to_detatch()
                survival_plot_spinner_locator = (
                    CohortComparisonLocators.LOADING_SPINNER_SURVIVAL_PLOT
                )
                self.wait_until_locator_is_detached(survival_plot_spinner_locator)
            except:
                self.wait_for_loading_spinner_to_detatch()

        if page_to_load == "set operations demo":
            # Need to wait for loading spinners to be present, for them to disappear
            try:
                self.wait_for_loading_spinner_to_be_visible(10000)
                self.wait_for_loading_spinner_to_detatch()
            except:
                self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
        self.wait_for_loading_spinner_table_to_detatch()
        self.wait_for_loading_spinner_to_detatch()
        self.wait_for_loading_spinner_cohort_bar_case_count_to_detatch()
