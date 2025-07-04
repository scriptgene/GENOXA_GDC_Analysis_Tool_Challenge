import React, { useEffect, useMemo, useRef, useState } from "react";
import { useDeepCompareEffect, useDeepCompareMemo } from "use-deep-compare";
import { ActionIcon, Drawer, ScrollArea } from "@mantine/core";
import { SortBy, useGetGenesQuery, useGetSsmsQuery } from "@gff/core";
import { humanify } from "src/utils";
import { SetData } from "./types";
import VerticalTable from "@/components/Table/VerticalTable";
import { createColumnHelper, SortingState } from "@tanstack/react-table";
import { LeftLongArrowIcon } from "@/utils/icons";

const PAGE_SIZE = 100;

interface SetDetailPanelProps {
  readonly set: SetData;
  readonly closePanel: () => void;
}

const SetDetailPanel: React.FC<SetDetailPanelProps> = ({
  set,
  closePanel,
}: SetDetailPanelProps) => {
  const [currentPage, setCurrentPage] = useState(0);
  const [tableData, setTableData] = useState([]);
  const [sortBy, setSortBy] = useState<SortBy[]>();
  const tableWrapperRef = useRef<HTMLDivElement>();
  const [scrollHeight, setScrollHeight] = useState(0);

  useEffect(() => {
    if (set === undefined) {
      setCurrentPage(0);
    }
  }, [set]);

  useEffect(() => {
    if (tableWrapperRef?.current?.clientHeight !== undefined) {
      setScrollHeight(tableWrapperRef?.current?.clientHeight);
    }
  }, [tableWrapperRef?.current?.clientHeight]);

  const {
    data: geneDetailData,
    isSuccess: isGeneSuccess,
    isFetching: isGeneFetching,
  } = useGetGenesQuery(
    {
      request: {
        filters: {
          op: "in",
          content: {
            field: "genes.gene_id",
            value: [`set_id:${set?.setId}`],
          },
        },
        fields: ["gene_id", "symbol"],
        size: PAGE_SIZE,
        from: currentPage * PAGE_SIZE,
        sortBy,
      },
      fetchAll: false,
    },
    { skip: set === undefined || set.setType === "ssms" },
  );

  const {
    data: ssmsDetailData,
    isSuccess: isMutationSuccess,
    isFetching: isMutationFetching,
  } = useGetSsmsQuery(
    {
      request: {
        filters: {
          op: "in",
          content: {
            field: "ssms.ssm_id",
            value: [`set_id:${set?.setId}`],
          },
        },
        fields: ["ssm_id"],
        expand: ["consequence.transcript", "consequence.transcript.gene"],
        size: PAGE_SIZE,
        from: currentPage * PAGE_SIZE,
        sortBy,
      },
      fetchAll: false,
    },
    { skip: set === undefined || set.setType === "genes" },
  );

  const responseData = useDeepCompareMemo(() => {
    if (set?.setType !== undefined) {
      if (set?.setType === "genes") {
        return isGeneSuccess && !isGeneFetching && geneDetailData?.hits
          ? [...geneDetailData.hits]
          : [];
      } else {
        return isMutationSuccess && !isMutationFetching
          ? ssmsDetailData?.hits?.map((ssm) => ({
              ssm_id: ssm.ssm_id,
              consequence: `${ssm?.consequence?.[0].transcript?.gene?.symbol} ${
                ssm?.consequence?.[0].transcript?.aa_change
                  ? ssm?.consequence?.[0].transcript?.aa_change
                  : ""
              } ${humanify({
                term: ssm?.consequence?.[0]?.transcript.consequence_type
                  ?.replace("_variant", "")
                  .replace("_", " "),
              })}`,
            }))
          : [];
      }
    }

    return [];
  }, [
    geneDetailData,
    isGeneFetching,
    isGeneSuccess,
    ssmsDetailData,
    isMutationFetching,
    isMutationSuccess,
    set?.setType,
  ]);

  // Append new data to existing table data
  useDeepCompareEffect(() => {
    setTableData([...tableData, ...responseData]);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [responseData]);

  useEffect(() => {
    if (set === undefined) {
      setTableData([]);
    }
  }, [set]);

  const setDetailPanelColumnHelper =
    createColumnHelper<(typeof tableData)[0]>();

  const setDetailPanelColumns = useMemo(
    () =>
      set?.setType === undefined
        ? []
        : set.setType === "genes"
        ? [
            setDetailPanelColumnHelper.accessor("gene_id", {
              id: "gene_id",
              header: "Gene ID",
            }),
            setDetailPanelColumnHelper.accessor("symbol", {
              id: "symbol",
              header: "Symbol",
            }),
          ]
        : [
            setDetailPanelColumnHelper.accessor("ssm_id", {
              id: "ssm_id",
              header: "Mutation ID",
            }),
            setDetailPanelColumnHelper.accessor("consequence", {
              id: "consequence",
              header: "Consequence",
              enableSorting: false,
            }),
          ],
    [set?.setType, setDetailPanelColumnHelper],
  );

  const scrollRef = useRef<HTMLDivElement>();

  const fetchMoreOnBottomReached = () => {
    if (scrollRef) {
      const { scrollHeight, scrollTop, clientHeight } = scrollRef.current;
      if (
        scrollHeight - scrollTop - clientHeight < 20 &&
        responseData.length > 0
      ) {
        setCurrentPage(currentPage + 1);
      }
    }
  };

  const displayData = useDeepCompareMemo(() => tableData, [tableData]);
  const [sorting, setSorting] = useState<SortingState>([]);

  useDeepCompareEffect(() => {
    setTableData([]);
    setCurrentPage(0);
    setSortBy(
      sorting.map((sort) => ({
        field: sort.id,
        direction: sort.desc ? "desc" : "asc",
      })),
    );
  }, [sorting]);

  return (
    <Drawer
      opened={set !== undefined}
      onClose={closePanel}
      position="right"
      classNames={{
        title: "w-full m-0",
        body: "h-full",
      }}
      title={
        <div
          data-testid="button-close-set-panel"
          className="flex flex-row gap-2 items-center w-full text-primary-darker font-bold p-2 border-b border-base-lighter"
        >
          <ActionIcon
            onClick={closePanel}
            aria-label="Close set panel"
            className="border-0"
          >
            <LeftLongArrowIcon size={30} className="text-primary-darker" />
          </ActionIcon>
          <>{set?.setName}</>
        </div>
      }
      size={"lg"}
      withCloseButton={false}
    >
      <div
        data-testid="table-set-information"
        className="h-full"
        ref={tableWrapperRef}
      >
        <ScrollArea
          h={scrollHeight}
          viewportRef={scrollRef}
          offsetScrollbars
          onScrollPositionChange={() => {
            fetchMoreOnBottomReached();
          }}
          className="pl-4"
        >
          <VerticalTable
            data={displayData}
            columns={setDetailPanelColumns}
            columnSorting="manual"
            sorting={sorting}
            status={
              isGeneSuccess || isMutationSuccess ? "fulfilled" : "pending"
            }
            setSorting={setSorting}
          />
        </ScrollArea>
      </div>
    </Drawer>
  );
};

export default SetDetailPanel;
