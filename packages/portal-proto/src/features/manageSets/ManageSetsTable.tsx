// This table can be found at /manage_sets
import React, { useMemo, useState } from "react";
import { useDeepCompareMemo, useDeepCompareEffect } from "use-deep-compare";
import { Checkbox, ActionIcon, Badge, Tooltip } from "@mantine/core";
import { showNotification } from "@mantine/notifications";
import { useCoreDispatch, removeSets, SetTypes } from "@gff/core";
import { createKeyboardAccessibleFunction } from "src/utils";
import download from "@/utils/download";
import useStandardPagination from "@/hooks/useStandardPagination";
import { SetData } from "./types";
import SetNameInput from "./SetNameInput";
import DeleteSetsNotification from "./DeleteSetsNotification";
import { HandleChangeInput } from "@/components/Table/types";
import VerticalTable from "@/components/Table/VerticalTable";
import { createColumnHelper, SortingState } from "@tanstack/react-table";
import { AlertIcon, DownloadIcon, TrashIcon } from "@/utils/icons";

interface CountBadgeProps {
  readonly count: number;
  readonly active: boolean;
  readonly openSetDetail: () => void;
}

const CountBadge: React.FC<CountBadgeProps> = ({
  count,
  active,
  openSetDetail,
}: CountBadgeProps) => {
  const disabled = count === 0;
  return (
    <Tooltip
      label="Set is either empty or deprecated."
      disabled={!disabled}
      withArrow
    >
      <Badge
        data-testid="text-set-count"
        variant={active ? "filled" : "outline"}
        radius="xs"
        className={`cursor-pointer w-16 ${
          active
            ? undefined
            : disabled
            ? "bg-nci-gray-lightest opacity-50 cursor-not-allowed"
            : "bg-white"
        }`}
        color={disabled ? "gray" : "primary"}
        leftSection={
          disabled ? <AlertIcon className="text-warningColor" /> : undefined
        }
        onClick={disabled ? undefined : openSetDetail}
        tabIndex={0}
        onKeyDown={
          disabled ? undefined : createKeyboardAccessibleFunction(openSetDetail)
        }
        aria-disabled={disabled}
        aria-label={`${count?.toLocaleString()} item${
          count > 1 ? "s" : ""
        } detailed list`}
        role="button"
      >
        {count?.toLocaleString() ?? "--"}
      </Badge>
    </Tooltip>
  );
};

interface ManageSetActionsProps {
  readonly set: SetData;
  readonly downloadType: string;
}

const ManageSetActions: React.FC<ManageSetActionsProps> = ({
  set,
  downloadType,
}: ManageSetActionsProps) => {
  const dispatch = useCoreDispatch();
  const { setId, setName, setType, count } = set;

  // add padding for both the icons
  return (
    <div className="flex flex-row items-center gap-3">
      <ActionIcon
        data-testid="button-delete-set"
        aria-label="Delete set"
        className="text-primary"
        onClick={() => {
          dispatch(removeSets([{ setId, setType }]));
          showNotification({
            message: <DeleteSetsNotification sets={[set]} />,
            closeButtonProps: { "aria-label": "Close notification" },
          });
        }}
        variant="outline"
      >
        <TrashIcon aria-hidden="true" />
      </ActionIcon>
      {count > 0 && (
        <ActionIcon
          data-testid="button-download-set"
          aria-label="Download set"
          className="text-primary"
          variant="outline"
          onClick={() => {
            download({
              endpoint: "tar_sets",
              params: {
                attachment: true,
                sets: [
                  {
                    id: setId,
                    type: downloadType,
                    filename: `${downloadType}_set_${setName.replace(
                      /[^A-Za-z0-9_.]/g,
                      "_",
                    )}.tsv`,
                  },
                ],
              },
              method: "POST",
              dispatch,
            });
          }}
        >
          <DownloadIcon aria-hidden="true" />
        </ActionIcon>
      )}
    </div>
  );
};

interface MangeSetsTableProps {
  readonly geneData: SetData[];
  readonly ssmData: SetData[];
  readonly setSelectedSets: (sets: SetData[]) => void;
  readonly detailSet: SetData;
  readonly setDetailSet: (set: SetData) => void;
}

type ManageSetsTableDataType = {
  setId: string;
  entityType: string;
  setName: string;
  count: number;
  set: SetData;
  setType: SetTypes;
};

const ManageSetsTable: React.FC<MangeSetsTableProps> = ({
  geneData,
  ssmData,
  setSelectedSets,
  detailSet,
  setDetailSet,
}) => {
  const tableData = useDeepCompareMemo(() => {
    return [
      ...ssmData.map((set) => {
        const { setName, count, setId, setType } = set;
        return {
          setId,
          entityType: "Mutations",
          setName,
          count,
          set,
          setType,
        };
      }),
      ...geneData.map((set) => {
        const { setName, count, setId, setType } = set;

        return {
          entityType: "Genes",
          setName,
          count,
          set,
          setId,
          setType,
        };
      }),
    ];
  }, [ssmData, geneData]);

  const getRowId = (originalRow: ManageSetsTableDataType) => {
    return originalRow.setId;
  };
  const [rowSelection, setRowSelection] = useState({});

  const pickedSets: SetData[] = Object.entries(rowSelection).reduce(
    (result, [id, isSelected]) => {
      if (isSelected) {
        const matchingItems = tableData.filter(
          (tableDatum) => tableDatum.setId === id,
        );
        result.push(...matchingItems);
      }
      return result;
    },
    [],
  );
  const [sorting, setSorting] = useState<SortingState>([]);

  useDeepCompareEffect(() => {
    setSelectedSets(pickedSets);
  }, [pickedSets, setSelectedSets]);

  const manageSetsTableColumnHelper =
    createColumnHelper<ManageSetsTableDataType>();

  const {
    handlePageChange,
    handlePageSizeChange,
    page,
    pages,
    size,
    from,
    total,
    displayedData,
  } = useStandardPagination(tableData);

  const handleChange = (obj: HandleChangeInput) => {
    switch (Object.keys(obj)?.[0]) {
      case "newPageSize":
        handlePageSizeChange(obj.newPageSize);
        break;
      case "newPageNumber":
        handlePageChange(obj.newPageNumber);
        break;
    }
  };

  const manageSetsColumn = useMemo(
    () => [
      manageSetsTableColumnHelper.display({
        id: "select",
        header: ({ table }) => (
          <Checkbox
            size="xs"
            classNames={{
              input: "checked:bg-accent checked:border-accent",
              label: "sr-only",
            }}
            {...{
              checked: table.getIsAllRowsSelected(),
              onChange: table.getToggleAllRowsSelectedHandler(),
            }}
            data-testid="checkbox-select-all-sets"
            label={`Select all sets on page ${page} of ${pages}`}
          />
        ),
        cell: ({ row }) => (
          <Checkbox
            size="xs"
            classNames={{
              input: "checked:bg-accent checked:border-accent",
            }}
            data-testid="checkbox-select-set"
            aria-label={row.original.setName}
            {...{
              checked: row.getIsSelected(),
              onChange: row.getToggleSelectedHandler(),
            }}
          />
        ),
        enableHiding: false,
      }),
      manageSetsTableColumnHelper.accessor("entityType", {
        id: "entityType",
        header: "Entity Type",
      }),
      manageSetsTableColumnHelper.accessor("setName", {
        id: "setName",
        header: "Name",
        cell: ({ getValue, row }) => (
          <SetNameInput
            setName={getValue()}
            setId={row.original.setId}
            setType={row.original.setType}
          />
        ),
      }),
      manageSetsTableColumnHelper.accessor("count", {
        id: "count",
        header: "# Items",
        cell: ({ getValue, row }) => (
          <CountBadge
            count={getValue()}
            active={detailSet?.setId === row.original.setId}
            openSetDetail={() => {
              setDetailSet({
                setId: row.original.setId,
                setName: row.original.setName,
                setType: row.original.setType,
                count: getValue(),
              });
            }}
          />
        ),
      }),
      manageSetsTableColumnHelper.display({
        id: "actions",
        header: "Actions",
        cell: ({ row }) => (
          <ManageSetActions
            set={row.original.set}
            downloadType={row.original.entityType === "Genes" ? "gene" : "ssm"}
          />
        ),
      }),
    ],
    [detailSet?.setId, setDetailSet, manageSetsTableColumnHelper, page, pages],
  );

  return (
    <div data-testid="table-manage-sets" className="w-full xl:w-3/4 pb-6">
      <VerticalTable
        data={displayedData}
        columns={manageSetsColumn}
        pagination={{
          page,
          pages,
          size,
          from,
          total,
          label: "set",
        }}
        handleChange={handleChange}
        columnSorting="enable"
        sorting={sorting}
        setSorting={setSorting}
        status="fulfilled"
        enableRowSelection={true}
        setRowSelection={setRowSelection}
        rowSelection={rowSelection}
        getRowId={getRowId}
      />
    </div>
  );
};

export default ManageSetsTable;
