import FunctionButton from "@/components/FunctionButton";
import { SaveCohortModal } from "@gff/portal-components";
import DarkFunctionButton from "@/components/StyledComponents/DarkFunctionButton";
import useStandardPagination from "@/hooks/useStandardPagination";
import {
  useCoreSelector,
  selectAvailableCohorts,
  FilterSet,
  selectCohortFilterSetById,
  fetchGdcCases,
  buildCohortGqlOperator,
  graphqlAPI,
} from "@gff/core";
import { Modal, Radio, Text, Loader } from "@mantine/core";
import { useMemo, useState } from "react";
import { MAX_CASE_IDS } from "./utils";
import { createColumnHelper } from "@tanstack/react-table";
import { HandleChangeInput } from "@/components/Table/types";
import VerticalTable from "@/components/Table/VerticalTable";
import { cohortActionsHooks } from "@/features/cohortBuilder/CohortManager/cohortActionHooks";
import { INVALID_COHORT_NAMES } from "@/features/cohortBuilder/utils";

export type WithOrWithoutCohortType = "with" | "without" | undefined;
export const SelectCohortsModal = ({
  opened,
  onClose,
  withOrWithoutCohort,
  pickedCases,
}: {
  opened: boolean;
  onClose: () => void;
  withOrWithoutCohort: WithOrWithoutCohortType;
  pickedCases: readonly string[];
}): JSX.Element => {
  const cohorts = useCoreSelector((state) => selectAvailableCohorts(state));
  const [checkedValue, setCheckedValue] = useState("");
  const cohortFilter = useCoreSelector((state) =>
    selectCohortFilterSetById(state, checkedValue),
  );
  const [saveCohortFilters, setSaveCohortFilters] = useState<
    FilterSet | undefined
  >(undefined);
  const [loading, setLoading] = useState(false);

  const isWithCohort = withOrWithoutCohort === "with";

  const cohortListData = useMemo(
    () =>
      cohorts
        ?.sort((a, b) => a.name.localeCompare(b.name))
        .map((cohort) => ({
          cohort_id: cohort.id,
          name: cohort.name,
          num_cases: cohort.counts.caseCount?.toLocaleString(),
        })),
    [cohorts],
  );

  const cohortListTableColumnHelper =
    createColumnHelper<(typeof cohortListData)[0]>();

  const cohortListTableColumn = useMemo(
    () => [
      cohortListTableColumnHelper.display({
        id: "select",
        header: "Select",
        cell: ({ row }) => (
          <Radio
            data-testid={`radio-${row.original.name}`}
            value={row.original.cohort_id}
            checked={checkedValue === row.original.cohort_id}
            onChange={(event) => {
              setCheckedValue(event.currentTarget.value);
            }}
          />
        ),
      }),
      cohortListTableColumnHelper.accessor("name", {
        id: "name",
        header: "Name",
        cell: ({ getValue, row }) => (
          <span data-testid={`text-${row.original.name}-cohort-name`}>
            {getValue()}
          </span>
        ),
      }),
      cohortListTableColumnHelper.accessor("num_cases", {
        id: "num_cases",
        header: "# Cases",
        cell: ({ getValue, row }) => (
          <span data-testid={`text-${row.original.name}-cohort-count`}>
            {getValue()}
          </span>
        ),
      }),
    ],
    [cohortListTableColumnHelper, checkedValue],
  );

  const {
    handlePageChange,
    handlePageSizeChange,
    page,
    pages,
    size,
    from,
    total,
    displayedData,
  } = useStandardPagination(cohortListData);

  const handleChange = (obj: HandleChangeInput) => {
    switch (Object.keys(obj)?.[0]) {
      case "newPageNumber":
        handlePageChange(obj.newPageNumber);
        break;
      case "newPageSize":
        handlePageSizeChange(obj.newPageSize);
        break;
    }
  };

  const createCohortFromCases = async () => {
    let resCases: string[];
    setLoading(true);

    try {
      const res = await fetchGdcCases({
        case_filters: buildCohortGqlOperator(cohortFilter),
        fields: ["case_id"],
        size: MAX_CASE_IDS,
      });
      resCases = res.data.hits.map((hit) => hit.case_id);
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
    } catch (error) {
      // TODO: how to handle this situation?
      // maybe show a modal and ask user to redo the task
    }

    const updatedCases = Array.from(
      new Set(
        !isWithCohort
          ? resCases.filter((id) => !pickedCases.includes(id))
          : pickedCases.concat(resCases),
      ),
    );

    const response = await graphqlAPI<any>(
      `mutation createSet(
      $input: CreateSetInput
    ) {
      sets {
        create {
          repository {
            case(input: $input) {
              set_id
              size
            }
          }
        }
      }
    }`,
      {
        input: {
          filters: {
            op: "and",
            content: [
              {
                op: "in",
                content: {
                  field: "cases.case_id",
                  value: updatedCases,
                },
              },
            ],
          },
          intent: "portal",
          set_type: "frozen",
        },
      },
    );
    const setId = response.data.sets.create.repository.case.set_id;

    setSaveCohortFilters({
      mode: "and",
      root: {
        "cases.case_id": {
          operator: "includes",
          field: "cases.case_id",
          operands: [`set_id:${setId}`],
        },
      },
    });
    setLoading(false);
  };

  const title = `save new cohort: existing cohort ${
    isWithCohort ? "with" : "without"
  } selected cases`;

  const description = `Select an existing cohort, then click Submit. This will save a new
    cohort that contains all the cases from your selected cohort ${
      isWithCohort ? "and" : "except"
    } the cases previously selected.`;

  const [showSaveCohort, setShowSaveCohorts] = useState(false);
  return (
    <>
      <Modal
        opened={opened}
        onClose={onClose}
        withCloseButton
        title={title}
        classNames={{
          content: "p-0 drop-shadow-lg",
          body: "flex flex-col justify-between min-h-[300px]",
        }}
        size="xl"
        zIndex={400}
      >
        <SaveCohortModal
          onClose={() => {
            setShowSaveCohorts(false);
            onClose();
          }}
          opened={showSaveCohort}
          filters={saveCohortFilters}
          hooks={cohortActionsHooks}
          invalidCohortNames={INVALID_COHORT_NAMES}
        />

        <div className="px-4">
          <Text className="text-sm mb-4 block font-content">{description}</Text>

          <VerticalTable
            customDataTestID="table-select-cohort"
            data={displayedData}
            columns={cohortListTableColumn}
            status="fulfilled"
            pagination={{
              page,
              pages,
              size,
              from,
              total,
              label: "cohort",
            }}
            handleChange={handleChange}
          />
        </div>
        <div
          data-testid="modal-button-container"
          className="bg-base-lightest flex p-4 gap-4 justify-end mt-4 rounded-b-lg sticky"
        >
          <FunctionButton data-testid="button-cancel" onClick={onClose}>
            Cancel
          </FunctionButton>
          <DarkFunctionButton
            data-testid="button-submit"
            disabled={!checkedValue}
            loading={loading}
            leftSection={
              loading ? <Loader size={15} color="white" /> : undefined
            }
            onClick={async () => {
              if (loading) return;
              await createCohortFromCases();
              setShowSaveCohorts(true);
            }}
          >
            Submit
          </DarkFunctionButton>
        </div>
      </Modal>
    </>
  );
};
