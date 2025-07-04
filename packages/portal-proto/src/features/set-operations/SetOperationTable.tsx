import VerticalTable from "@/components/Table/VerticalTable";
import {
  ColumnDef,
  createColumnHelper,
  SortingState,
} from "@tanstack/react-table";
import { useState, useId } from "react";
import { useDeepCompareMemo } from "use-deep-compare";
import CountButtonWrapperForSetsAndCases from "./CountButtonWrapperForSetsAndCases";
import { Checkbox, Tooltip } from "@mantine/core";
import { createSetFiltersByKey, ENTITY_TYPE_TO_CREATE_SET_HOOK } from "./utils";
import { pickBy } from "lodash";
import { SelectedEntities } from "./types";
import DownloadButtonTotal from "./DownloadButton";
import { GqlUnion, useSetOperationsCasesTotalQuery } from "@gff/core";

type SetOperationTableDataType = {
  setOperation: string;
  count: number;
  operationKey: string;
};

const setOperationTableColumnsHelper =
  createColumnHelper<SetOperationTableDataType>();

export const SetOperationTable = ({
  data,
  selectedSets,
  setSelectedSets,
  queryHook,
  entityType,
  sets,
}: {
  readonly sets: SelectedEntities;
  readonly entityType: "cohort" | "genes" | "mutations";
  readonly data: {
    readonly label: string;
    readonly key: string;
    readonly value: number;
  }[];
  readonly queryHook: typeof useSetOperationsCasesTotalQuery;
  selectedSets: {
    [k: string]: boolean;
  };
  setSelectedSets: React.Dispatch<
    React.SetStateAction<{
      [k: string]: boolean;
    }>
  >;
}): JSX.Element => {
  const unionFilter = {
    op: "or",
    content: Object.keys(pickBy(selectedSets, (v) => v)).map((set) =>
      createSetFiltersByKey(set, entityType, sets),
    ),
  } as GqlUnion;

  const { data: totalSelectedSets, isFetching } = queryHook({
    filters: {
      filters: unionFilter,
    },
  });
  const totalCount =
    Object.keys(pickBy(selectedSets, (v) => v)).length > 0
      ? totalSelectedSets
      : 0;
  const [rowSelection, setRowSelection] = useState({});
  const [operationTableSorting, setOperationTableSorting] =
    useState<SortingState>([]);
  const componentId = useId();

  const setOperationTableData: SetOperationTableDataType[] = useDeepCompareMemo(
    () =>
      data.map((r) => ({
        setOperation: r.label,
        count: r.value,
        operationKey: r.key,
      })),
    [data],
  );

  const setOperationTableColumns = useDeepCompareMemo<
    ColumnDef<SetOperationTableDataType>[]
  >(
    () => [
      {
        id: "select",
        header: "Select",
        cell: ({ row }) => (
          <Tooltip
            label="This region contains 0 items"
            disabled={row.original.count !== 0}
            position="right"
          >
            <Checkbox
              data-testid={`checkbox-${row.original.operationKey}-set-operations`}
              size="xs"
              classNames={{
                input: "checked:bg-accent checked:border-accent",
              }}
              value={row.original.operationKey}
              id={`${componentId}-setOperation-${row.original.operationKey}`}
              checked={selectedSets[row.original.operationKey]}
              onChange={(e) => {
                setSelectedSets({
                  ...selectedSets,
                  [e.target.value]: !selectedSets[e.target.value],
                });
              }}
              disabled={row.original.count === 0}
            />
          </Tooltip>
        ),
      },
      setOperationTableColumnsHelper.accessor("setOperation", {
        header: "Set Operation",
        enableSorting: false,
        cell: ({ getValue, row }) => (
          <label
            data-testid={`text-${row.original.operationKey}-label-set-operations`}
            htmlFor={`${componentId}-setOperation-${row.original.operationKey}`}
          >
            {getValue()}
          </label>
        ),
      }),
      setOperationTableColumnsHelper.accessor("count", {
        header: "# Items",
        cell: ({ row, getValue }) => (
          <CountButtonWrapperForSetsAndCases
            count={getValue()}
            filters={createSetFiltersByKey(
              row.original.operationKey,
              entityType,
              sets,
            )}
            entityType={entityType}
          />
        ),
        enableSorting: true,
      }),
      setOperationTableColumnsHelper.display({
        id: "download",
        header: "Download",
        cell: ({ row }) => (
          <DownloadButtonTotal
            filters={createSetFiltersByKey(
              row.original.operationKey,
              entityType,
              sets,
            )}
            entityType={entityType}
            setKey={row.original.setOperation}
            createSetHook={ENTITY_TYPE_TO_CREATE_SET_HOOK[entityType]}
            disabled={row.original.count === 0}
          />
        ),
      }),
    ],

    [selectedSets, entityType, sets, setSelectedSets, componentId],
  );

  return (
    <VerticalTable
      customDataTestID="table-set-operations"
      data={setOperationTableData}
      columns={setOperationTableColumns}
      enableRowSelection={true}
      setRowSelection={setRowSelection}
      rowSelection={rowSelection}
      showControls={false}
      status="fulfilled"
      sorting={operationTableSorting}
      setSorting={setOperationTableSorting}
      columnSorting="enable"
      customAriaLabel="Overlap Table"
      footer={
        <tr data-testid="row-union-of-selected-sets">
          <td />
          <td className="p-2 font-bold">Union of selected sets:</td>
          <td className="w-52 px-2.5">
            <CountButtonWrapperForSetsAndCases
              count={isFetching ? 0 : totalCount}
              filters={unionFilter}
              entityType={entityType}
            />
          </td>
          <td className="p-2.5">
            <DownloadButtonTotal
              setKey="union-of"
              filters={unionFilter}
              createSetHook={ENTITY_TYPE_TO_CREATE_SET_HOOK[entityType]}
              entityType={entityType}
              disabled={totalCount === 0}
            />
          </td>
        </tr>
      }
    />
  );
};
