import React from "react";
import { useDeepCompareMemo } from "use-deep-compare";
import Link from "next/link";
import { Loader } from "@mantine/core";
import { useGetAnnotationsQuery, useQuickSearchQuery } from "@gff/core";
import { SummaryHeader } from "@/components/Summary/SummaryHeader";
import { HeaderTitle } from "@/components/tailwindComponents";
import { HorizontalTable } from "@/components/HorizontalTable";
import { SummaryErrorHeader } from "@/components/Summary/SummaryErrorHeader";
import { PencilSquareIcon } from "@/utils/icons";

interface AnnotationSummaryProps {
  readonly annotationId: string;
}

const AnnotationSummary: React.FC<AnnotationSummaryProps> = ({
  annotationId,
}) => {
  const {
    data: annotationData,
    isSuccess,
    isFetching,
  } = useGetAnnotationsQuery({
    request: {
      filters: {
        op: "=",
        content: {
          field: "annotation_id",
          value: annotationId,
        },
      },
      expand: ["project"],
    },
  });

  const annotation = annotationData?.hits?.[0];

  const { data: entityData } = useQuickSearchQuery(annotation?.entity_id);

  const entityLink = useDeepCompareMemo(() => {
    if (annotation === undefined) {
      return undefined;
    }

    if (
      entityData?.searchList === undefined ||
      entityData?.searchList.length == 0
    ) {
      return annotation?.entity_id ?? "--";
    } else {
      if (annotation?.entity_type === "case") {
        return (
          <Link
            data-testid="link-entity-annotation-summary"
            href={`/cases/${annotation?.entity_id}`}
            className="underline text-utility-link"
          >
            {annotation?.entity_id}
          </Link>
        );
      } else {
        const entityType = atob(entityData?.searchList?.[0].id).split(":")[0];
        if (entityType === "File") {
          return (
            <Link
              data-testid="link-entity-annotation-summary"
              href={`/files/${annotation?.entity_id}`}
              className="underline text-utility-link"
            >
              {annotation?.entity_id}
            </Link>
          );
        } else {
          return (
            <Link
              data-testid="link-entity-annotation-summary"
              href={`/cases/${annotation?.case_id}?bioId=${annotation?.entity_id}`}
              className="underline text-utility-link"
            >
              {annotation?.entity_id}
            </Link>
          );
        }
      }
    }
  }, [entityData, annotation]);

  const idTableData = useDeepCompareMemo(
    () => [
      { headerName: "Annotation UUID", values: [annotationId] },
      { headerName: "Entity UUID", values: [entityLink] },
      {
        headerName: "Entity ID",
        values: [annotation?.entity_submitter_id ?? "--"],
      },
      { headerName: "Entity Type", values: [annotation?.entity_type ?? "--"] },
      {
        headerName: "Case UUID",
        values: [
          annotation?.case_id ? (
            annotation.classification === "Redaction" &&
            annotation.entity_type === "case" ? (
              annotation.case_id
            ) : (
              <Link
                data-testid="link-case-annotation-summary"
                href={`/cases/${annotation.case_id}`}
                className="underline text-utility-link"
                key={`case_link_${annotation.case_id}`}
              >
                {annotation.case_id}
              </Link>
            )
          ) : (
            "--"
          ),
        ],
      },
      {
        headerName: "Case ID",
        values: [annotation?.case_submitter_id ?? "--"],
      },
    ],
    [annotationId, annotation, entityLink],
  );

  const tableData = useDeepCompareMemo(
    () => [
      {
        headerName: "Project",
        values: [
          annotation?.project?.project_id ? (
            <Link
              data-testid="link-project-annotation-summary"
              href={`/projects/${annotation?.project?.project_id}`}
              className="underline text-utility-link"
              key={`project_link_${annotation?.project?.project_id}`}
            >
              {annotation?.project?.project_id}
            </Link>
          ) : (
            "--"
          ),
        ],
      },
      {
        headerName: "Classification",
        values: [annotation?.classification ?? "--"],
      },
      { headerName: "Category", values: [annotation?.category ?? "--"] },
      {
        headerName: "Created On",
        values: [annotation?.created_datetime ?? "--"],
      },
      { headerName: "Status", values: [annotation?.status ?? "--"] },
    ],
    [annotation],
  );

  return isFetching ? (
    <Loader />
  ) : isSuccess && annotationData?.hits.length === 0 ? (
    <SummaryErrorHeader label="Annotation Not Found" />
  ) : (
    <>
      <SummaryHeader
        Icon={PencilSquareIcon}
        headerTitleLeft="Annotation"
        headerTitle={annotationId}
      />
      <div className="flex flex-col gap-2 mx-4 mt-6 mb-4">
        <HeaderTitle>Summary</HeaderTitle>
        <div
          data-testid="table-summary-annotation-summary"
          className="flex mb-8"
        >
          <div className="basis-1/2">
            <HorizontalTable
              customDataTestID="table-left-summary-annotation-summary"
              tableData={idTableData}
            />
          </div>
          <div className="basis-1/2">
            <HorizontalTable
              customDataTestID="table-right-summary-annotation-summary"
              tableData={tableData}
            />
          </div>
        </div>
        <HeaderTitle>Notes</HeaderTitle>
        <p
          data-testid="table-notes-annotation-summary"
          className="border-1 border-base-lighter bg-primary-content-lightest p-2 font-content"
        >
          {annotation?.notes}
        </p>
      </div>
    </>
  );
};

export default AnnotationSummary;
