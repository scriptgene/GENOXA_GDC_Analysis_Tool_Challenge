import { Button } from "@mantine/core";
import { saveAs } from "file-saver";
import {
  FilterSet,
  buildCohortGqlOperator,
  useCoreSelector,
  selectCurrentCohortFilters,
} from "@gff/core";
import { DropdownWithIcon } from "@/components/DropdownWithIcon/DropdownWithIcon";
import { CasesCohortButtonFromFilters } from "@/features/cases/CasesView/CasesCohortButton";
import { useIsDemoApp } from "@/hooks/useIsDemoApp";
import { getFormattedTimestamp } from "@/utils/date";
import {
  CategoricalBins,
  CustomInterval,
  DisplayData,
  NamedFromTo,
  SelectedFacet,
} from "../types";
import { DEMO_COHORT_FILTERS } from "../constants";
import {
  createFiltersFromSelectedValues,
  formatPercent,
  useDataDimension,
} from "../utils";
import { useMemo } from "react";
import { DropdownIcon } from "@/utils/icons";

interface CardControlsProps {
  readonly continuous: boolean;
  readonly field: string;
  readonly fieldName: string;
  readonly displayedData: DisplayData;
  readonly yTotal: number;
  readonly setBinningModalOpen: (open: boolean) => void;
  readonly customBinnedData: CategoricalBins | NamedFromTo[] | CustomInterval;
  readonly setCustomBinnedData:
    | ((bins: CategoricalBins) => void)
    | ((bins: NamedFromTo[] | CustomInterval) => void);
  readonly selectedFacets: SelectedFacet[];
  readonly dataDimension?: string;
}

const CardControls: React.FC<CardControlsProps> = ({
  continuous,
  field,
  fieldName,
  displayedData,
  yTotal,
  setBinningModalOpen,
  customBinnedData,
  setCustomBinnedData,
  selectedFacets,
  dataDimension,
}: CardControlsProps) => {
  const isDemoMode = useIsDemoApp();
  const displayDataDimension = useDataDimension(field);

  const downloadTSVFile = () => {
    const header = [
      displayDataDimension ? `${fieldName} (${dataDimension})` : fieldName,
      "# Cases",
    ];
    const body = displayedData.map(({ displayName, count }) =>
      [displayName, `${count} (${formatPercent(count, yTotal)})`].join("\t"),
    );
    const tsv = [header.join("\t"), body.join("\n")].join("\n");

    saveAs(
      new Blob([tsv], {
        type: "text/tsv",
      }),
      `${field.split(".").at(-1)}-table.${getFormattedTimestamp()}.tsv`,
    );
  };

  const cohortFilters = useCoreSelector((state) =>
    isDemoMode ? DEMO_COHORT_FILTERS : selectCurrentCohortFilters(state),
  );

  const filters: FilterSet = useMemo(
    () =>
      createFiltersFromSelectedValues(
        continuous,
        field,
        selectedFacets,
        customBinnedData,
      ),
    [continuous, field, selectedFacets, customBinnedData],
  );

  return (
    <>
      <div className="flex justify-between gap-2 py-2 flex-reverse-wrap md:flex-wrap">
        <div className="flex flex-wrap gap-2">
          <div className="order-2 md:order-1">
            <CasesCohortButtonFromFilters
              filters={
                selectedFacets.length === 0
                  ? undefined
                  : buildCohortGqlOperator(filters)
              }
              case_filters={
                selectedFacets.length === 0
                  ? undefined
                  : buildCohortGqlOperator(cohortFilters)
              }
              numCases={
                selectedFacets.length === 0
                  ? 0
                  : selectedFacets
                      .map((facet) => facet.numCases)
                      .reduce((a, b) => a + b)
              }
            />
          </div>
          <Button
            data-testid="button-tsv-cdave-card"
            className="bg-base-max text-primary border-primary order-1 md:order-2"
            onClick={downloadTSVFile}
          >
            TSV
          </Button>
        </div>
        <div className="flex items-end">
          <DropdownWithIcon
            customTargetButtonDataTestId="button-customize-bins"
            RightSection={
              <div className="border-l pl-1 -mr-2">
                <DropdownIcon size={20} aria-hidden="true" />
              </div>
            }
            TargetButtonChildren={"Customize Bins"}
            disableTargetWidth={true}
            dropdownElements={[
              {
                title: "Edit Bins",
                onClick: () => setBinningModalOpen(true),
              },
              {
                title: "Reset to Default",
                disabled: customBinnedData === null,
                onClick: () => setCustomBinnedData(null),
              },
            ]}
          />
        </div>
      </div>
    </>
  );
};

export default CardControls;
