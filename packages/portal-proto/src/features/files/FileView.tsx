import React, { useRef, useState } from "react";
import { useDeepCompareMemo } from "use-deep-compare";
import {
  GdcFile,
  useCoreSelector,
  selectCart,
  Modals,
  selectCurrentModal,
} from "@gff/core";
import { get } from "lodash";
import dynamic from "next/dynamic";
import fileSize from "filesize";
import tw from "tailwind-styled-components";
import { AddToCartButton, RemoveFromCartButton } from "../cart/updateCart";
import {
  formatDataForHorizontalTable,
  mapGdcFileToCartFile,
  parseSlideDetailsInfo,
  shouldDisplayReferenceGenome,
} from "./utils";
import { BAMSlicingModal } from "@/components/Modals/BAMSlicingModal/BAMSlicingModal";
import { NoAccessToProjectModal } from "@/components/Modals/NoAccessToProjectModal";
import { BAMSlicingButton } from "@/features/files/BAMSlicingButton";
import { DownloadFile } from "@/components/DownloadButtons";
import { AgreementModal } from "@/components/Modals/AgreementModal";
import { fileInCart } from "src/utils";
import { GeneralErrorModal } from "@/components/Modals/GeneraErrorModal";
import { HeaderTitle } from "@/components/tailwindComponents";
import { SummaryCard } from "@/components/Summary/SummaryCard";
import { SummaryHeader } from "@/components/Summary/SummaryHeader";
import GenericLink from "@/components/GenericLink";
import AssociatedCB from "./AssociatedCB";
import DownstreamAnalyses from "./DownstreamAnalyses";
import SourceFiles from "./SourceFiles";
import ReadGroups from "./ReadGroups";
import FileVersions from "./FileVersions";
import AnnnotationsTable from "./AnnotationsTable";
import FilesIcon from "public/user-flow/icons/summary/files.svg";
import { useSynchronizedRowHeights } from "@/components/HorizontalTable/useSynchronizedRowHeights";

interface LeftSideElementForHeaderProps {
  readonly isFileInCart: boolean;
  readonly file: GdcFile;
  readonly bamActive: boolean;
  readonly setFileToDownload: React.Dispatch<React.SetStateAction<GdcFile>>;
}

const LeftSideElementForHeader: React.FC<LeftSideElementForHeaderProps> = ({
  isFileInCart,
  file,
  bamActive,
  setFileToDownload,
}: LeftSideElementForHeaderProps) => (
  <div className="flex gap-4">
    {!isFileInCart ? (
      <AddToCartButton files={mapGdcFileToCartFile([file])} />
    ) : (
      <RemoveFromCartButton files={mapGdcFileToCartFile([file])} />
    )}
    {file.data_format === "BAM" &&
      file.data_type === "Aligned Reads" &&
      file?.index_files?.length > 0 && (
        <BAMSlicingButton isActive={bamActive} file={file} />
      )}

    <DownloadFile
      customDataTestID="button-download-file-summary"
      inactiveText="Download"
      activeText="Processing"
      file={file}
      displayVariant="header-subtle"
      setfileToDownload={setFileToDownload}
    />
  </div>
);

const ImageViewer = dynamic(
  () => import("../../components/ImageViewer/ImageViewer"),
  {
    ssr: false,
  },
);

export interface FileViewProps {
  readonly file?: GdcFile;
  readonly isModal?: boolean;
}

const DivWithMargin = tw.div`mt-8 flex flex-col gap-2`;

export const FileView: React.FC<FileViewProps> = ({
  file,
  isModal,
}: FileViewProps) => {
  const currentCart = useCoreSelector((state) => selectCart(state));
  const modal = useCoreSelector((state) => selectCurrentModal(state));
  const [bamActive, setBamActive] = useState(false);
  const [fileToDownload, setFileToDownload] = useState(file);
  const isFileInCart = fileInCart(currentCart, file.file_id);
  const shouldDisplayRefGenome = shouldDisplayReferenceGenome(file);
  const leftAnalysisTableRef = useRef<HTMLTableElement>(null);
  const rightAnalysisTableRef = useRef<HTMLTableElement>(null);
  const leftBAMMetricsTableRef = useRef<HTMLTableElement>(null);
  const rightBAMMetricsTableRef = useRef<HTMLTableElement>(null);
  useSynchronizedRowHeights([leftAnalysisTableRef, rightAnalysisTableRef]);
  useSynchronizedRowHeights([leftBAMMetricsTableRef, rightBAMMetricsTableRef]);

  const formattedDataForFileProperties = useDeepCompareMemo(
    () =>
      formatDataForHorizontalTable(file, [
        {
          field: "file_name",
          name: "Name",
          modifier: (v) => <span className="break-all">{v}</span>,
        },
        {
          field: "access",
          name: "Access",
        },
        {
          field: "id",
          name: "UUID",
        },
        {
          field: "data_format",
          name: "Data Format",
        },
        {
          field: "file_size",
          name: "Size",
          modifier: fileSize,
        },
        {
          field: "md5sum",
          name: "MD5 Checksum",
        },
        {
          field: "project_id",
          name: "Project",
          modifier: (v) => <GenericLink path={`/projects/${v}`} text={v} />,
        },
      ]),
    [file],
  );

  const formattedDataForDataInformation = useDeepCompareMemo(
    () =>
      formatDataForHorizontalTable(file, [
        {
          field: "data_category",
          name: "Data Category",
        },
        {
          field: "data_type",
          name: "Data Type",
        },
        {
          field: "experimental_strategy",
          name: "Experimental Strategy",
        },
        {
          field: "platform",
          name: "Platform",
        },
      ]),
    [file],
  );

  const formattedDataForAnalysis = useDeepCompareMemo(
    () =>
      formatDataForHorizontalTable(file, [
        {
          field: "analysis.workflow_type",
          name: "Workflow Type",
        },
      ]),
    [file],
  );

  const formattedDataForQCMetrics = useDeepCompareMemo(
    () =>
      formatDataForHorizontalTable(file, [
        {
          field: "total_reads",
          name: "Total Reads",
        },
        {
          field: "average_insert_size",
          name: "Average Insert Size",
        },
        {
          field: "average_read_length",
          name: "Average Read Length",
        },
        {
          field: "average_base_quality",
          name: "Average Base Quality",
        },
        {
          field: "mean_coverage",
          name: "Mean Coverage",
        },
        {
          field: "pairs_on_diff_chr",
          name: "Pairs On Diff Chr",
        },
        {
          field: "contamination",
          name: "Contamination",
        },
        {
          field: "contamination_error",
          name: "Contamination Error",
        },

        {
          field: "proportion_reads_mapped",
          name: "Proportion Reads Mapped",
        },
        {
          field: "proportion_reads_duplicated",
          name: "Proportion Reads Duplicated",
        },
        {
          field: "proportion_base_mismatch",
          name: "Proportion Base Mismatch",
        },
        {
          field: "proportion_targets_no_coverage",
          name: "Proportion Targets No Coverage",
        },
        {
          field: "proportion_coverage_10x",
          name: "Proportion Coverage 10X",
        },
        {
          field: "proportion_coverage_30x",
          name: "Proportion Coverage 30X",
        },
        {
          field: "msi_score",
          name: "MSI Score",
        },
        {
          field: "msi_status",
          name: "MSI Status",
        },
      ]),
    [file],
  );

  const isAllMetricsValUndefined = formattedDataForQCMetrics.every(
    (metric) => metric.values.length === 1 && metric.values.includes("--"),
  );

  return (
    <div className="relative">
      <SummaryHeader
        Icon={FilesIcon}
        headerTitleLeft="File"
        headerTitle={file.file_name}
        isModal={isModal}
        leftElement={
          <LeftSideElementForHeader
            file={file}
            isFileInCart={isFileInCart}
            bamActive={bamActive}
            setFileToDownload={setFileToDownload}
          />
        }
      />
      <div className={`${!isModal ? "mt-6" : "mt-4"} mx-4`}>
        <div className="flex flex-col lg:flex-row gap-8 xl:gap-4">
          <div className="flex-1">
            <SummaryCard
              customDataTestID="table-file-properties-file-summary"
              title="File Properties"
              tableData={formattedDataForFileProperties}
            />
          </div>
          <div className="flex-1">
            <SummaryCard
              customDataTestID="table-data-information-file-summary"
              tableData={formattedDataForDataInformation}
              title="Data Information"
            />
          </div>
        </div>

        {get(file, "data_type") === "Slide Image" && (
          <DivWithMargin data-testid="table-slide-image-file-summary">
            <HeaderTitle>Slide Image Viewer</HeaderTitle>
            <ImageViewer
              imageId={file?.file_id}
              tableData={parseSlideDetailsInfo(file)}
            />
          </DivWithMargin>
        )}
        <DivWithMargin data-testid="table-associated-cases-biospecimens-file-summary">
          {file?.associated_entities?.length > 0 ? (
            <AssociatedCB
              cases={file?.cases}
              associated_entities={file?.associated_entities}
            />
          ) : (
            <>
              <HeaderTitle>Associated Cases/Biospecimens</HeaderTitle>
              <h3 className="p-2 mx-4 text-primary-content-darker">
                No cases or biospecimen found.
              </h3>
            </>
          )}
        </DivWithMargin>

        {file.data_format === "BAM" &&
          file.data_type === "Aligned Reads" &&
          !isAllMetricsValUndefined && (
            <div className="flex mt-8">
              <div className="basis-1/2">
                <SummaryCard
                  customDataTestID="table-left-bam-metrics-file-summary"
                  title="BAM Metrics"
                  tableData={formattedDataForQCMetrics.slice(0, 8)}
                  ref={leftBAMMetricsTableRef}
                  enableSync={true}
                />
              </div>

              <div className="basis-1/2">
                <SummaryCard
                  customDataTestID="table-right-bam-metrics-file-summary"
                  title="" // should be empty
                  tableData={formattedDataForQCMetrics.slice(8, 16)}
                  ref={rightBAMMetricsTableRef}
                  enableSync={true}
                />
              </div>
            </div>
          )}

        {file?.analysis && (
          <>
            <div className="mt-8 flex flex-col lg:flex-row gap-8 xl:gap-4">
              <div
                className={`flex-1 ${!shouldDisplayRefGenome ? "flex" : ""}`}
              >
                <div className="basis-1/2">
                  <SummaryCard
                    customDataTestID="table-analysis-file-summary"
                    title="Analysis"
                    tableData={formattedDataForAnalysis}
                    ref={leftAnalysisTableRef}
                    enableSync={!shouldDisplayRefGenome}
                  />
                </div>
              </div>

              {shouldDisplayRefGenome && (
                <div className="flex-1">
                  <SummaryCard
                    customDataTestID="table-reference-genome-file-summary"
                    title="Reference Genome"
                    tableData={[
                      { headerName: "Genome Build", values: ["GRCh38.p0"] },
                      { headerName: "Genome Name", values: ["GRCh38.d1.vd1"] },
                    ]}
                  />
                </div>
              )}
            </div>
            {file?.analysis?.input_files?.length > 0 && (
              <DivWithMargin>
                <SourceFiles
                  inputFiles={file.analysis.input_files}
                  currentCart={currentCart}
                  setFileToDownload={setFileToDownload}
                />
              </DivWithMargin>
            )}

            {file?.analysis?.metadata && (
              <DivWithMargin data-testid="table-read-groups-file-summary">
                <ReadGroups readGroups={file?.analysis.metadata.read_groups} />
              </DivWithMargin>
            )}
          </>
        )}
        {file?.downstream_analyses?.some(
          (byWorkflowType) => byWorkflowType?.output_files?.length > 0,
        ) && (
          <DivWithMargin data-testid="table-downstream-analyses-files-file-summary">
            <DownstreamAnalyses
              downstream_analyses={file?.downstream_analyses}
              currentCart={currentCart}
              setFileToDownload={setFileToDownload}
            />
          </DivWithMargin>
        )}
        <FileVersions file_id={file.file_id} />
        {file?.annotations?.length > 0 && (
          <div
            className={`mb-16 ${isModal ? "scroll-mt-36" : "scroll-mt-72"}`}
            id="annotations"
          >
            <AnnnotationsTable annotations={file.annotations} />
          </div>
        )}

        {modal === Modals.NoAccessToProjectModal && (
          <NoAccessToProjectModal openModal />
        )}
        {modal === Modals.BAMSlicingModal && (
          <BAMSlicingModal openModal file={file} setActive={setBamActive} />
        )}

        {modal === Modals.GeneralErrorModal && <GeneralErrorModal openModal />}

        {modal === Modals.AgreementModal && (
          <AgreementModal
            openModal
            file={fileToDownload}
            dbGapList={fileToDownload.acl}
          />
        )}
      </div>
    </div>
  );
};
