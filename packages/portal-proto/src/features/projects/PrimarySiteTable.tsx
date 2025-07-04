import React, { useState, useEffect, useMemo } from "react";
import { useGetProjectsPrimarySitesAllQuery } from "@gff/core";
import useStandardPagination from "@/hooks/useStandardPagination";
import CohortCreationButton from "@/components/CohortCreationButton";
import { useDeepCompareMemo } from "use-deep-compare";
import {
  ColumnOrderState,
  ExpandedState,
  Row,
  SortingState,
  VisibilityState,
  createColumnHelper,
} from "@tanstack/react-table";
import VerticalTable from "@/components/Table/VerticalTable";
import { HandleChangeInput } from "@/components/Table/types";
import ExpandRowComponent from "@/components/Table/ExpandRowComponent";
import { HeaderTitle } from "@/components/tailwindComponents";
import SubrowPrimarySiteDiseaseType from "@/components/SubrowPrimarySiteDiseaseType/SubrowPrimarySiteDiseaseType";
import { ArraySeparatedSpan } from "@/components/ArraySeparatedSpan/ArraySeparatedSpan";
import TotalItems from "@/components/Table/TotalItem";

interface PrimarySiteTableProps {
  readonly projectId: string;
  readonly primarySites: string[];
}

const PrimarySiteTable: React.FC<PrimarySiteTableProps> = ({
  projectId,
  primarySites,
}: PrimarySiteTableProps) => {
  const [filteredTableData, setFilteredTableData] = useState([]);
  const [searchTerm, setSearchTerm] = useState<string>("");

  const { data, isFetching } = useGetProjectsPrimarySitesAllQuery({
    projectId,
    primary_sites: primarySites,
  });

  const formattedData = useDeepCompareMemo(
    () =>
      !isFetching
        ? data?.map((datum) => ({
            primary_site: datum.primary_site,
            disease_type: datum.disease_types,
            cases: datum.casesTotal,
            experimental_strategy: datum.files__experimental_strategy,
            files: datum.filesTotal?.toLocaleString(),
          }))
        : [],
    [isFetching, data],
  );

  const primarySitesTableColumnHelper = useMemo(
    () => createColumnHelper<(typeof formattedData)[0]>(),
    [],
  );

  const primarySitesTableColumns = useMemo(
    () => [
      primarySitesTableColumnHelper.accessor("primary_site", {
        id: "primary_site",
        header: "Primary Site",
      }),
      primarySitesTableColumnHelper.accessor("disease_type", {
        id: "disease_type",
        header: "Disease Type",
        cell: ({ row, getValue }) => (
          <ExpandRowComponent
            value={getValue()}
            title="Disease Types"
            isRowExpanded={row.getIsExpanded()}
          />
        ),
      }),
      primarySitesTableColumnHelper.accessor("cases", {
        id: "cases",
        header: "Cases",
        cell: ({ row }) => (
          <CohortCreationButton
            label={row.original.cases?.toLocaleString()}
            numCases={row.original.cases}
            filters={{
              mode: "and",
              root: {
                "cases.project.project_id": {
                  field: "cases.project.project_id",
                  operator: "includes",
                  operands: [projectId],
                },
                "cases.primary_site": {
                  field: "cases.primary_site",
                  operator: "includes",
                  operands: [row.original.primary_site.toLowerCase()],
                },
              },
            }}
          />
        ),
        meta: {
          sortingFn: (rowA, rowB) => {
            if (rowA.cases > rowB.cases) {
              return 1;
            }
            if (rowA.cases < rowB.cases) {
              return -1;
            }
            return 0;
          },
        },
      }),
      primarySitesTableColumnHelper.accessor("experimental_strategy", {
        id: "experimental_strategy",
        header: "Experimental Strategy",
        cell: ({ getValue }) => (
          <ArraySeparatedSpan data={getValue() as string[]} />
        ),
        enableSorting: false,
      }),
      primarySitesTableColumnHelper.accessor("files", {
        id: "files",
        header: "Files",
        meta: {
          sortingFn: (rowA, rowB) => {
            const a = Number(rowA.files.replace(/,/g, ""));
            const b = Number(rowB.files.replace(/,/g, ""));
            if (a > b) return 1;
            if (a < b) return -1;
            return 0;
          },
        },
      }),
    ],
    [primarySitesTableColumnHelper, projectId],
  );

  useEffect(() => {
    if (searchTerm) {
      setFilteredTableData(
        formattedData.filter((obj) => {
          return obj.primary_site?.toLowerCase().indexOf(searchTerm) > -1;
        }),
      );
    } else {
      setFilteredTableData(formattedData);
    }
  }, [searchTerm, formattedData]);

  const [columnOrder, setColumnOrder] = useState<ColumnOrderState>(
    primarySitesTableColumns.map((column) => column.id as string), //must start out with populated columnOrder so we can splice
  );
  const [columnVisibility, setColumnVisibility] = useState<VisibilityState>({});
  const [expanded, setExpanded] = useState<ExpandedState>({});
  const getRowId = (originalRow: (typeof formattedData)[0]) => {
    return originalRow.primary_site;
  };
  const [expandedRowId, setExpandedRowId] = useState(null);
  const [sorting, setSorting] = useState<SortingState>([
    {
      id: "primary_site",
      desc: false,
    },
  ]);

  const {
    handlePageChange,
    handlePageSizeChange,
    handleSortByChange,
    displayedData,
    ...paginationProps
  } = useStandardPagination(filteredTableData, primarySitesTableColumns);

  useEffect(() => handleSortByChange(sorting), [sorting, handleSortByChange]);

  const handleChange = (obj: HandleChangeInput) => {
    switch (Object.keys(obj)?.[0]) {
      case "newPageSize":
        handlePageSizeChange(obj.newPageSize);
        break;
      case "newPageNumber":
        handlePageChange(obj.newPageNumber);
        break;
      case "newSearch":
        setSearchTerm(obj.newSearch?.toLowerCase());
        break;
    }
  };

  const handleExpand = (
    row: Row<(typeof formattedData)[0]>,
    columnId: string,
  ) => {
    if (
      Object.keys(expanded).length > 0 &&
      row.original.primary_site === expandedRowId
    ) {
      setExpanded({});
    } else if ((row.original[columnId] as string[]).length > 1) {
      setExpanded({ [row.original.primary_site]: true });
      setExpandedRowId(row.original.primary_site);
    }
  };

  return (
    <VerticalTable
      customDataTestID="table-primary-sites-project-summary"
      tableTotalDetail={
        <TotalItems total={paginationProps?.total} itemName="primary site" />
      }
      data={displayedData}
      columns={primarySitesTableColumns}
      additionalControls={<HeaderTitle>Primary Sites</HeaderTitle>}
      columnSorting="manual"
      sorting={sorting}
      setSorting={setSorting}
      expanded={expanded}
      getRowCanExpand={() => true}
      setExpanded={handleExpand}
      expandableColumnIds={["disease_type"]}
      renderSubComponent={({ row, clickedColumnId }) => (
        <SubrowPrimarySiteDiseaseType row={row} columnId={clickedColumnId} />
      )}
      setColumnVisibility={setColumnVisibility}
      columnVisibility={columnVisibility}
      columnOrder={columnOrder}
      setColumnOrder={setColumnOrder}
      pagination={{
        ...paginationProps,
        label: "Primary Site",
      }}
      status={isFetching ? "pending" : "fulfilled"}
      search={{
        enabled: true,
        tooltip: "e.g. Bronchus and lung, Kidney",
      }}
      showControls={true}
      handleChange={handleChange}
      baseZIndex={300}
      getRowId={getRowId}
    />
  );
};

export default PrimarySiteTable;
