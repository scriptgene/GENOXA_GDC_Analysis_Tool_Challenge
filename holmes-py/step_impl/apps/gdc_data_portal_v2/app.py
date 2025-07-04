import os

from step_impl.base.base_page import BasePage as Shared
from ..gdc_data_portal_v2.pages.header_section import HeaderSection
from step_impl.apps.gdc_data_portal_v2.pages.home_page import HomePage
from step_impl.apps.gdc_data_portal_v2.pages.analysis_center_page import (
    AnalysisCenterPage,
)
from step_impl.apps.gdc_data_portal_v2.pages.clinical_data_analysis import (
    ClinicalDataAnalysisPage,
)
from step_impl.apps.gdc_data_portal_v2.pages.modal import Modal
from ..gdc_data_portal_v2.pages.repository_page import RepositoryPage
from ..gdc_data_portal_v2.pages.cohort_builder_page import CohortBuilderPage
from ..gdc_data_portal_v2.pages.file_summary_page import FileSummaryPage
from ..gdc_data_portal_v2.pages.case_summary_page import CaseSummaryPage
from ..gdc_data_portal_v2.pages.cohort_bar import CohortBar
from ..gdc_data_portal_v2.pages.projects_page import ProjectsPage
from ..gdc_data_portal_v2.pages.mutation_frequency_page import MutationFrequencyPage
from ..gdc_data_portal_v2.pages.manage_sets_page import ManageSetsPage
from ..gdc_data_portal_v2.pages.cohort_comparison_page import CohortComparisonPage
from ..gdc_data_portal_v2.pages.project_summary_page import ProjectSummaryPage
from ..gdc_data_portal_v2.pages.set_operations_page import SetOperationsPage
from ..gdc_data_portal_v2.pages.cohort_case_view import CohortCaseViewPage
from ..gdc_data_portal_v2.pages.cart_page import CartPage
from ..gdc_data_portal_v2.pages.gene_summary_page import GeneSummaryPage
from ..gdc_data_portal_v2.pages.mutation_summary_page import MutationSummaryPage
from ..gdc_data_portal_v2.pages.browse_annotations import BrowseAnnotationPage


class GDCDataPortalV2App:
    def __init__(self, webdriver):  # webdriver is page now.
        app_endpoint_var = (
            "APP_ENDPOINT_PROD"
            if not os.getenv("APP_ENVIRONMENT")
            else f"APP_ENDPOINT_{os.environ['APP_ENVIRONMENT']}"
        )
        self.URL = f"{os.getenv(app_endpoint_var)}"
        self.driver = webdriver
        self.init_pages()

    def navigate(self):
        self.driver.goto(self.URL)

    def init_pages(self):
        # 'Shared' contains common functions and locators seen throughout the Data Portal.
        # It uses the code contained in base_page.py
        self.shared = Shared(self.driver)
        self.header_section = HeaderSection(self.driver, self.URL)
        self.modal = Modal(self.driver, self.URL)
        self.home_page = HomePage(self.driver, self.URL)
        self.repository_page = RepositoryPage(self.driver, self.URL)
        self.cohort_builder_page = CohortBuilderPage(self.driver, self.URL)
        self.analysis_center_page = AnalysisCenterPage(self.driver, self.URL)
        self.clinical_data_analysis = ClinicalDataAnalysisPage(self.driver, self.URL)
        self.file_summary_page = FileSummaryPage(self.driver, self.URL)
        self.case_summary_page = CaseSummaryPage(self.driver, self.URL)
        self.cohort_bar = CohortBar(self.driver, self.URL)
        self.projects_page = ProjectsPage(self.driver, self.URL)
        self.mutation_frequency_page = MutationFrequencyPage(self.driver, self.URL)
        self.manage_sets_page = ManageSetsPage(self.driver, self.URL)
        self.cohort_comparison_page = CohortComparisonPage(self.driver, self.URL)
        self.project_summary_page = ProjectSummaryPage(self.driver, self.URL)
        self.set_operations_page = SetOperationsPage(self.driver, self.URL)
        self.cohort_case_view_page = CohortCaseViewPage(self.driver, self.URL)
        self.cart_page = CartPage(self.driver, self.URL)
        self.gene_summary_page = GeneSummaryPage(self.driver, self.URL)
        self.mutation_summary_page = MutationSummaryPage(self.driver, self.URL)
        self.browse_annotations = BrowseAnnotationPage(self.driver, self.URL)
