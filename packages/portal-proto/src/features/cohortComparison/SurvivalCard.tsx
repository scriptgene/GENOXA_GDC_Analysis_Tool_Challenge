import { Alert, LoadingOverlay, Paper, Tooltip } from "@mantine/core";
import {
  buildCohortGqlOperator,
  useGetSurvivalPlotQuery,
  FilterSet,
  GqlIntersection,
  GqlOperation,
  useCreateCaseSetFromFiltersMutation,
} from "@gff/core";
import SurvivalPlot from "../charts/SurvivalPlot/SurvivalPlot";
import makeIntersectionFilters from "./makeIntersectionFilters";
import CohortCreationButton from "@/components/CohortCreationButton";
import { SurvivalPlotTypes } from "../charts/SurvivalPlot/types";
import { useDeepCompareEffect } from "use-deep-compare";
import { emptySurvivalPlot } from "../genomic/types";

const survivalDataCompletenessFilters: GqlOperation[] = [
  {
    op: "or",
    content: [
      {
        op: "and",
        content: [
          {
            op: ">",
            content: {
              field: "demographic.days_to_death",
              value: 0,
            },
          },
        ],
      },
      {
        op: "and",
        content: [
          {
            op: ">",
            content: {
              field: "diagnoses.days_to_last_follow_up",
              value: 0,
            },
          },
        ],
      },
    ],
  },
  {
    op: "not",
    content: { field: "demographic.vital_status" },
  },
];

export const makeSurvivalCaseFilters = (
  primaryCohortSetId: string,
  comparisonCohortSetId: string,
): GqlIntersection => ({
  op: "and",
  content: [
    ...survivalDataCompletenessFilters,
    {
      op: "and",
      content: [
        {
          op: "in",
          content: {
            field: "cases.case_id",
            value: [`set_id:${primaryCohortSetId}`],
          },
        },
        {
          op: "excludeifany",
          content: {
            field: "cases.case_id",
            value: `set_id:${comparisonCohortSetId}`,
          },
        },
      ],
    },
  ],
});

const tooltipLabel = (
  <>
    <p>
      Criteria for including a case from your cohorts in the survival analysis:
    </p>
    <p>- the case is in only one of your cohorts, not both</p>
    <p>- the case has the required data for the analysis</p>
  </>
);

interface SurvivalCardProps {
  readonly counts: number[];
  readonly cohorts?: {
    primary_cohort: {
      filter: FilterSet;
      name: string;
    };
    comparison_cohort: {
      filter: FilterSet;
      name: string;
    };
  };
  readonly setSurvivalPlotSelectable: (selectable: boolean) => void;
  readonly caseSetIds: string[];
  readonly isSetsloading: boolean;
}

const SurvivalCard: React.FC<SurvivalCardProps> = ({
  counts,
  cohorts,
  setSurvivalPlotSelectable,
  caseSetIds,
  isSetsloading,
}: SurvivalCardProps) => {
  const filters = makeIntersectionFilters(
    buildCohortGqlOperator(cohorts?.primary_cohort.filter),
    buildCohortGqlOperator(cohorts?.comparison_cohort.filter),
    caseSetIds,
  );
  const [createSet] = useCreateCaseSetFromFiltersMutation();

  const generatePrimaryFilters = async () => {
    return await createSet({
      filters: makeSurvivalCaseFilters(caseSetIds[0], caseSetIds[1]),
      intent: "portal",
      set_type: "frozen",
    })
      .unwrap()
      .then((setId) => {
        return {
          mode: "and",
          root: {
            "cases.case_id": {
              field: "cases.case_id",
              operands: [`set_id:${setId}`],
              operator: "includes",
            },
          },
        } as FilterSet;
      });
  };

  const generateComparisonFilters = async () => {
    return await createSet({
      filters: makeSurvivalCaseFilters(caseSetIds[1], caseSetIds[0]),
      intent: "portal",
      set_type: "frozen",
    })
      .unwrap()
      .then((setId) => {
        return {
          mode: "and",
          root: {
            "cases.case_id": {
              field: "cases.case_id",
              operands: [`set_id:${setId}`],
              operator: "includes",
            },
          },
        } as FilterSet;
      });
  };

  const { data, isUninitialized, isFetching, isError } =
    useGetSurvivalPlotQuery({
      filters: [filters.cohort1, filters.cohort2],
    });

  useDeepCompareEffect(() => {
    setSurvivalPlotSelectable(data?.survivalData.length !== 0);
  }, [data, setSurvivalPlotSelectable]);

  const cohort1Count = data?.survivalData[0]
    ? data.survivalData[0].donors?.length
    : 0;
  const cohort2Count = data?.survivalData[1]
    ? data.survivalData[1].donors?.length
    : 0;

  const isLoading = isSetsloading || isFetching || isUninitialized;

  return (
    <Paper data-testid="card-analysis-survival-cohort-comparison" p="md">
      <h2 className="font-heading text-lg font-semibold">Survival Analysis</h2>
      {data?.survivalData.length === 0 ? (
        <div className="p-1">
          No Survival data available for this Cohort Comparison
        </div>
      ) : (
        <>
          {isError ? (
            <Alert>Something`&apos;s gone wrong</Alert>
          ) : (
            <div className="relative">
              <LoadingOverlay
                visible={isLoading}
                zIndex={0}
                data-testid="loading-spinner"
              />
              <SurvivalPlot
                plotType={SurvivalPlotTypes.cohortComparison}
                data={isLoading ? emptySurvivalPlot : data}
                hideLegend
                noDataMessage="No Survival data available for this Cohort Comparison"
                isLoading={isLoading}
              />
            </div>
          )}
          <div className="font-heading mt-[1.5rem]">
            <table className="bg-base-max w-full text-left text-base-contrast-max border-base-light border-1">
              <thead>
                <tr className="bg-base-lightest border-b-base-light border-b-2">
                  <th>
                    <Tooltip label={tooltipLabel}>
                      <span className="underline decoration-dashed pl-2">
                        Cases included in Analysis
                      </span>
                    </Tooltip>
                  </th>
                  <th>
                    # Cases S<sub>1</sub>
                  </th>
                  <th>
                    % Cases S<sub>1</sub>
                  </th>
                  <th>
                    # Cases S<sub>2</sub>
                  </th>
                  <th>
                    % Cases S<sub>2</sub>
                  </th>
                </tr>
              </thead>
              <tbody
                data-testid="text-analysis-overall survival analysis"
                className="font-content text-md"
              >
                <tr>
                  <td className="pl-2">Overall Survival Analysis</td>
                  <td>
                    {isLoading ? (
                      "..."
                    ) : (
                      <CohortCreationButton
                        numCases={cohort1Count}
                        label={cohort1Count.toLocaleString()}
                        filtersCallback={generatePrimaryFilters}
                      />
                    )}
                  </td>
                  <td>
                    {isLoading
                      ? "..."
                      : `${((cohort1Count / counts[0]) * 100).toFixed(0)}%`}
                  </td>
                  <td>
                    {isLoading ? (
                      "..."
                    ) : (
                      <CohortCreationButton
                        numCases={cohort2Count}
                        label={cohort2Count.toLocaleString()}
                        filtersCallback={generateComparisonFilters}
                      />
                    )}
                  </td>
                  <td>
                    {isLoading
                      ? "..."
                      : `${((cohort2Count / counts[1]) * 100).toFixed(0)}%`}
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
        </>
      )}
    </Paper>
  );
};

export default SurvivalCard;
