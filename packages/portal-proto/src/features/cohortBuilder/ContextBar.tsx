import React, { useState } from "react";
import { useDeepCompareEffect } from "use-deep-compare";
import { CollapsibleContainer } from "@/components/CollapsibleContainer";
import { Loader, Tabs } from "@mantine/core";
import { ContextualCasesView } from "../cases/CasesView/CasesView";
import CohortCountButton from "./CohortCountButton";
import CohortManager from "./CohortManager/CohortManager";
import {
  useCoreDispatch,
  useCoreSelector,
  selectAvailableCohorts,
  selectCurrentCohortId,
  setActiveCohort,
  Modals,
  showModal,
  selectCurrentCohortFilters,
  buildCohortGqlOperator,
  useGetAllFilesMutation,
  GdcFile,
  selectCart,
  GqlOperation,
  useCurrentCohortCounts,
} from "@gff/core";
import download from "src/utils/download";
import SummaryFacets from "./SummaryFacets";
import { SecondaryTabStyle } from "@/features/cohortBuilder/style";
import { mapGdcFileToCartFile } from "@/features/files/utils";
import { addToCart } from "@/features/cart/updateCart";
import { DropdownWithIcon } from "@/components/DropdownWithIcon/DropdownWithIcon";
import {
  INITIAL_SUMMARY_FIELDS,
  METADATA_EXPAND_PROPS,
  METADATA_FIELDS,
  SAMPLE_SHEET_FIELDS,
} from "./utils";
import StickyControl from "./StickyControl";
import { useViewportSize } from "@mantine/hooks";
import {
  SummaryChartIcon,
  TableIcon,
  CohortFilterIcon,
  DownloadIcon,
} from "@/utils/icons";
import { getFormattedTimestamp } from "@/utils/date";

const ContextBar = ({
  handleIsSticky,
  isSticky,
  isFetchingCohorts,
}: {
  handleIsSticky: (isSticky: boolean) => void;
  isSticky: boolean;
  isFetchingCohorts: boolean;
}): JSX.Element => {
  const coreDispatch = useCoreDispatch();
  const cohorts = useCoreSelector((state) => selectAvailableCohorts(state));
  const [summaryFields] = useState(INITIAL_SUMMARY_FIELDS);
  const [activeTab, setActiveTab] = useState<string | null>("summary");
  const [isGroupCollapsed, setIsGroupCollapsed] = useState(true);
  const { width } = useViewportSize();

  /* download active */
  const [downloadManifestActive, setDownloadManifestActive] = useState(false);
  const [downloadMetadataActive, setDownloadMetadataActive] = useState(false);
  const [downloadSampleSheetActive, setDownloadSampleSheetActive] =
    useState(false);
  const [biospecimenDownloadActive, setBiospecimenDownloadActive] =
    useState(false);
  const [clinicalDownloadActive, setClinicalDownloadActive] = useState(false);
  /* download active end */

  const currentCohortId = useCoreSelector((state) =>
    selectCurrentCohortId(state),
  );

  useDeepCompareEffect(() => {
    if (currentCohortId === undefined && cohorts.length > 0) {
      coreDispatch(setActiveCohort(cohorts[0].id));
    }
  }, [currentCohortId, cohorts, coreDispatch]);

  const cohortFilters = useCoreSelector((state) =>
    selectCurrentCohortFilters(state),
  );
  const cohortCounts = useCurrentCohortCounts();
  const caseCounts = cohortCounts?.data?.caseCount;

  const currentCart = useCoreSelector((state) => selectCart(state));
  const [getFileSizeSliceData, { isLoading: fetchingFilesForCohort }] =
    useGetAllFilesMutation();

  const getAllCohortFiles = (callback, filters) => {
    getFileSizeSliceData({ caseFilters: filters })
      .unwrap()
      .then((data: GdcFile[]) => {
        return mapGdcFileToCartFile(data);
      })
      .then((cartFiles) => {
        callback(cartFiles, currentCart, coreDispatch);
      });
  };

  const downloadFilter: GqlOperation =
    buildCohortGqlOperator(cohortFilters) ?? ({} as GqlOperation);

  const handleMetadataDownload = () => {
    setDownloadMetadataActive(true);
    download({
      endpoint: "files",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        case_filters: downloadFilter,
        fields: METADATA_FIELDS.join(","),
        format: "JSON",
        pretty: "True",
        attachment: "True",
        expand: METADATA_EXPAND_PROPS.join(","),
        filename: `metadata.cohort.${new Date()
          .toISOString()
          .slice(0, 10)}.json`,
      },
      done: () => setDownloadMetadataActive(false),
    });
  };

  const handleManifestDownload = () => {
    setDownloadManifestActive(true);
    download({
      endpoint: "files",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        case_filters: downloadFilter,
        return_type: "manifest",
        size: 10000,
        filename: `gdc_manifest.${getFormattedTimestamp({
          includeTimes: true,
        })}.txt`,
      },
      done: () => setDownloadManifestActive(false),
    });
  };

  const handleSampleSheetDownload = () => {
    setDownloadSampleSheetActive(true);
    download({
      endpoint: "files",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        case_filters: downloadFilter,
        tsv_format: "sample-sheet",
        fields: SAMPLE_SHEET_FIELDS.join(","),
        format: "tsv",
        pretty: "True",
        attachment: "True",
        filename: `gdc_sample_sheet.${new Date()
          .toISOString()
          .slice(0, 10)}.tsv`,
      },
      done: () => setDownloadSampleSheetActive(false),
    });
  };

  const handleClinicalTSVDownload = () => {
    setClinicalDownloadActive(true);
    download({
      endpoint: "clinical_tar",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        filename: `clinical.cohort.${new Date()
          .toISOString()
          .slice(0, 10)}.tar.gz`,
        case_filters: downloadFilter,
        size: caseCounts,
      },
      done: () => setClinicalDownloadActive(false),
    });
  };

  const handleClinicalJSONDownload = () => {
    setClinicalDownloadActive(true);
    download({
      endpoint: "clinical_tar",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        format: "JSON",
        pretty: true,
        filename: `clinical.cohort.${new Date()
          .toISOString()
          .slice(0, 10)}.json`,
        filters: downloadFilter,
        size: caseCounts,
      },
      done: () => setClinicalDownloadActive(false),
    });
  };

  const handleBiospeciemenTSVDownload = () => {
    setBiospecimenDownloadActive(true);
    download({
      endpoint: "biospecimen_tar",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        filename: `biospecimen.cohort.${new Date()
          .toISOString()
          .slice(0, 10)}.tar.gz`,
        filters: downloadFilter,
        size: caseCounts,
      },
      done: () => setBiospecimenDownloadActive(false),
    });
  };

  const handleBiospeciemenJSONDownload = () => {
    setBiospecimenDownloadActive(true);
    download({
      endpoint: "biospecimen_tar",
      method: "POST",
      dispatch: coreDispatch,
      params: {
        format: "JSON",
        pretty: true,
        filename: `biospecimen.cohort.${new Date()
          .toISOString()
          .slice(0, 10)}.json`,
        filters: downloadFilter,
        size: caseCounts,
      },
      done: () => setBiospecimenDownloadActive(false),
    });
  };

  return (
    <CollapsibleContainer
      Top={<CohortManager isFetchingCohorts={isFetchingCohorts} />}
      isCollapsed={isGroupCollapsed}
      toggle={() => setIsGroupCollapsed(!isGroupCollapsed)}
      onlyIcon={false}
      isContextBar={true}
      tooltipText={
        isGroupCollapsed
          ? "Access additional cases details and features"
          : "Collapse additional cases details and features"
      }
      tooltipPosition="left"
      TargetElement={
        <CohortCountButton
          customDataTestID="button-cases-cohort-bar"
          countName="caseCount"
          singularLabel="case"
          bold
          capitalize
        />
      }
      ExtraControl={
        <StickyControl isSticky={isSticky} handleIsSticky={handleIsSticky} />
      }
    >
      <div
        className={`bg-nci-violet-lightest ${isGroupCollapsed && "invisible"}`}
        aria-hidden={isGroupCollapsed}
      >
        <div className="relative p-4">
          <div className="flex gap-1 pb-4 relative lg:absolute lg:pb-0">
            <DropdownWithIcon
              customTargetButtonDataTestId="button-files-cases-summary"
              dropdownElements={[
                {
                  title: "Add to Cart",
                  onClick: () =>
                    getAllCohortFiles(
                      addToCart,
                      buildCohortGqlOperator(cohortFilters),
                    ),
                  icon: fetchingFilesForCohort ? (
                    <Loader size={14} />
                  ) : undefined,
                },
                {
                  title: "Download Manifest",
                  onClick: handleManifestDownload,
                  icon: downloadManifestActive ? (
                    <Loader size={14} />
                  ) : undefined,
                },
                {
                  title: "Metadata",
                  onClick: handleMetadataDownload,
                  icon: downloadMetadataActive ? (
                    <Loader size={14} />
                  ) : undefined,
                },
                {
                  title: "Sample Sheet",
                  onClick: handleSampleSheetDownload,
                  icon: downloadSampleSheetActive ? (
                    <Loader size={14} />
                  ) : undefined,
                },
              ]}
              TargetButtonChildren={
                <CohortCountButton countName="fileCount" singularLabel="file" />
              }
              targetButtonTooltip={
                cohortCounts?.data?.fileCount === 0
                  ? "No files in current cohort"
                  : undefined
              }
              LeftSection={
                <DownloadIcon
                  size="1rem"
                  aria-hidden="true"
                  className="hidden md:block"
                />
              }
              targetButtonDisabled={
                cohortCounts.status !== "fulfilled" ||
                cohortCounts?.data?.fileCount === 0
              }
            />

            <DropdownWithIcon
              customTargetButtonDataTestId="button-custom-filters-cases-summary"
              dropdownElements={[
                {
                  title: "Cases",
                  onClick: () =>
                    coreDispatch(
                      showModal({ modal: Modals.GlobalCaseSetModal }),
                    ),
                },
                {
                  title: "Mutated Genes",
                  onClick: () =>
                    coreDispatch(
                      showModal({ modal: Modals.GlobalGeneSetModal }),
                    ),
                },
                {
                  title: "Somatic Mutations",
                  onClick: () =>
                    coreDispatch(
                      showModal({ modal: Modals.GlobalMutationSetModal }),
                    ),
                },
              ]}
              TargetButtonChildren="Custom Filters"
              LeftSection={
                <CohortFilterIcon
                  size="1rem"
                  aria-hidden="true"
                  className="hidden md:block"
                />
              }
              menuLabelText="Filter your cohort by:"
              menuLabelCustomClass="font-bold text-primary"
              targetButtonDisabled={cohortCounts.status !== "fulfilled"}
            />

            {activeTab === "summary" && (
              <>
                <DropdownWithIcon
                  customTargetButtonDataTestId="button-biospecimen-cases-summary"
                  dropdownElements={[
                    {
                      title: "JSON ",
                      onClick: handleBiospeciemenJSONDownload,
                      icon: <DownloadIcon aria-label="Download" />,
                    },
                    {
                      title: "TSV",
                      onClick: handleBiospeciemenTSVDownload,
                      icon: <DownloadIcon aria-label="Download" />,
                    },
                  ]}
                  TargetButtonChildren={
                    biospecimenDownloadActive ? "Processing" : "Biospecimen"
                  }
                  LeftSection={
                    <span className="hidden md:block">
                      {biospecimenDownloadActive ? (
                        <Loader size={20} />
                      ) : (
                        <DownloadIcon size="1rem" aria-hidden="true" />
                      )}
                    </span>
                  }
                  targetButtonDisabled={cohortCounts.status !== "fulfilled"}
                />

                <DropdownWithIcon
                  customTargetButtonDataTestId="button-clinical-cases-summary"
                  dropdownElements={[
                    {
                      title: "JSON",
                      onClick: handleClinicalJSONDownload,
                      icon: <DownloadIcon aria-label="Download" />,
                    },
                    {
                      title: "TSV",
                      onClick: handleClinicalTSVDownload,
                      icon: <DownloadIcon aria-label="Download" />,
                    },
                  ]}
                  TargetButtonChildren={
                    clinicalDownloadActive ? "Processing" : "Clinical"
                  }
                  LeftSection={
                    <span className="hidden md:block">
                      {biospecimenDownloadActive ? (
                        <Loader size={20} />
                      ) : (
                        <DownloadIcon size="1rem" aria-hidden="true" />
                      )}
                    </span>
                  }
                  targetButtonDisabled={cohortCounts.status !== "fulfilled"}
                />
              </>
            )}
          </div>
          <Tabs
            data-testid="filters-facets-summary-view"
            variant="pills"
            classNames={{
              tab: SecondaryTabStyle,
              list: "mb-4 gap-0",
              root: "border-0",
              panel: "h-max",
            }}
            data-tour="cohort_summary"
            defaultValue="summary"
            keepMounted={false}
            value={activeTab}
            onChange={setActiveTab}
          >
            <Tabs.List justify={width < 1024 ? "flex-start" : "flex-end"}>
              <Tabs.Tab
                data-testid="button-summary-view"
                data-tour="cohort_summary_charts"
                value="summary"
                leftSection={<SummaryChartIcon aria-hidden="true" />}
              >
                Summary View
              </Tabs.Tab>

              <Tabs.Tab
                data-testid="button-table-view"
                data-tour="cohort_summary_table"
                value="table"
                leftSection={<TableIcon aria-hidden="true" />}
              >
                Table View
              </Tabs.Tab>
            </Tabs.List>
            <Tabs.Panel value="summary">
              {!isGroupCollapsed && ( //dont load unless shown TODO address this in CollapsibleContainer
                <SummaryFacets fields={summaryFields} />
              )}
            </Tabs.Panel>
            <Tabs.Panel value="table">
              <ContextualCasesView />
            </Tabs.Panel>
          </Tabs>
        </div>
      </div>
    </CollapsibleContainer>
  );
};

export default ContextBar;
