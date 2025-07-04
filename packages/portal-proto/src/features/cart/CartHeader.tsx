import React, { useState } from "react";
import {
  CartFile,
  CartSummaryData,
  showModal,
  useCoreDispatch,
  CoreDispatch,
  useCoreSelector,
  selectCurrentModal,
  Modals,
  useFetchUserDetailsQuery,
} from "@gff/core";
import fileSize from "filesize";
import { Button, Loader, Menu } from "@mantine/core";
import CartSizeLimitModal from "@/components/Modals/CartSizeLimitModal";
import CartDownloadModal from "@/components/Modals/CartDownloadModal";
import { DownloadButton } from "@/components/DownloadButtons";
import download from "src/utils/download";
import { removeFromCart } from "./updateCart";
import { focusStyles } from "@/utils/index";
import { cartAboveLimit } from "./utils";
import {
  ArrowDropDownIcon,
  CartIcon,
  DownloadIcon,
  FileIcon,
  PersonIcon,
  SaveIcon,
} from "@/utils/icons";
import { getFormattedTimestamp } from "@/utils/date";

const buttonStyle =
  "bg-base-max text-primary border-primary data-disabled:opacity-50 data-disabled:bg-base-max data-disabled:text-primary";

const downloadCart = (
  filesByCanAccess: Record<string, CartFile[]>,
  dbGapList: string[],
  setActive: (active: boolean) => void,
  dispatch: CoreDispatch,
) => {
  if (
    cartAboveLimit(filesByCanAccess) ||
    (filesByCanAccess?.false || []).length > 0 ||
    dbGapList.length > 0
  ) {
    dispatch(showModal({ modal: Modals.CartDownloadModal }));
    setActive(false);
  } else {
    download({
      endpoint: "data",
      method: "POST",
      dispatch,
      params: {
        ids: filesByCanAccess.true.map((file) => file.file_id),
        annotations: true,
        related_files: true,
      },
      done: () => setActive(false),
    });
  }
};

const downloadManifest = (
  cart: CartFile[],
  setActive: (active: boolean) => void,
  dispatch: CoreDispatch,
) => {
  download({
    endpoint: "files",
    method: "POST",
    dispatch,
    params: {
      filters: {
        op: "in",
        content: {
          field: "files.file_id",
          value: cart.map((file) => file.file_id),
        },
      },
      return_type: "manifest",
      size: 10000,
      filename: `gdc_manifest.${getFormattedTimestamp({
        includeTimes: true,
      })}.txt`,
    },
    done: () => setActive(false),
  });
};

interface CartHeaderProps {
  summaryData: CartSummaryData;
  cart: CartFile[];
  filesByCanAccess: Record<string, CartFile[]>;
  dbGapList: string[];
}

const CartHeader: React.FC<CartHeaderProps> = ({
  summaryData,
  cart,
  filesByCanAccess,
  dbGapList,
}: CartHeaderProps) => {
  const dispatch = useCoreDispatch();
  const { data: userDetails } = useFetchUserDetailsQuery();
  const [cartDownloadActive, setCartDownloadActive] = useState(false);
  const [manifestDownloadActive, setManifestDownloadActive] = useState(false);
  const [clinicalTSVDownloadActive, setClinicalTSVDownloadActive] =
    useState(false);
  const [clinicalJSONDownloadActive, setClinicalJSONDownloadActive] =
    useState(false);
  const [biospecimenTSVDownloadActive, setBiospecimenTSVDownloadActive] =
    useState(false);
  const [biospecimenJSONDownloadActive, setBiospecimenJSONDownloadActive] =
    useState(false);
  const [metadataDownloadActive, setMetadataDownloadActive] = useState(false);
  const [sampleSheetDownloadActive, setSampleSheetDownloadActive] =
    useState(false);
  const modal = useCoreSelector((state) => selectCurrentModal(state));

  return (
    <>
      {modal === Modals.CartSizeLimitModal && <CartSizeLimitModal openModal />}
      {modal === Modals.CartDownloadModal && (
        <CartDownloadModal
          openModal
          user={userDetails?.data}
          filesByCanAccess={filesByCanAccess}
          dbGapList={dbGapList}
          setActive={setCartDownloadActive}
        />
      )}
      <div
        className="bg-primary text-primary-contrast-darkest flex flex-col-reverse 2xl:flex-row 2xl:items-center gap-4 w-full p-4"
        data-testid="cart-header"
      >
        <div className="flex flex-wrap gap-2">
          <Menu width="target">
            <Menu.Target>
              <Button
                data-testid="button-download-cart"
                classNames={{
                  root: `${buttonStyle} ${focusStyles}`,
                }}
                leftSection={
                  <DownloadIcon
                    aria-hidden="true"
                    size="1rem"
                    className="hidden xl:block"
                  />
                }
                rightSection={
                  <div className="border-l pl-1 -mr-2">
                    <ArrowDropDownIcon size="1.5em" aria-hidden="true" />
                  </div>
                }
              >
                Download Cart
              </Button>
            </Menu.Target>
            <Menu.Dropdown data-testid="dropdown-menu-options">
              <Menu.Item
                onClick={() => {
                  setManifestDownloadActive(true);
                  downloadManifest(cart, setManifestDownloadActive, dispatch);
                }}
                leftSection={
                  manifestDownloadActive ? (
                    <Loader size={15} />
                  ) : (
                    <DownloadIcon aria-hidden="true" />
                  )
                }
              >
                Manifest
              </Menu.Item>
              <Menu.Item
                onClick={() => {
                  setCartDownloadActive(true);
                  downloadCart(
                    filesByCanAccess,
                    dbGapList,
                    setCartDownloadActive,
                    dispatch,
                  );
                }}
                leftSection={
                  cartDownloadActive ? (
                    <Loader size={15} />
                  ) : (
                    <DownloadIcon aria-hidden="true" />
                  )
                }
              >
                Cart
              </Menu.Item>
            </Menu.Dropdown>
          </Menu>
          <Menu width="target">
            <Menu.Target>
              <Button
                data-testid="button-download-associated-data"
                classNames={{
                  root: `${buttonStyle} ${focusStyles}`,
                }}
                leftSection={
                  <DownloadIcon
                    aria-hidden="true"
                    size="1rem"
                    className="hidden xl:block"
                  />
                }
                rightSection={
                  <div className="border-l pl-1 -mr-2">
                    <ArrowDropDownIcon size="1.5em" aria-hidden="true" />
                  </div>
                }
              >
                Biospecimen
              </Button>
            </Menu.Target>
            <Menu.Dropdown data-testid="dropdown-menu-options">
              <Menu.Item
                component={DownloadButton}
                classNames={{ item: "font-normal" }}
                displayVariant="subtle"
                activeText="Processing"
                inactiveText="JSON"
                preventClickEvent
                showIcon={true}
                endpoint="biospecimen_tar"
                setActive={setBiospecimenJSONDownloadActive}
                active={biospecimenJSONDownloadActive}
                filename={`biospecimen.cart.${new Date()
                  .toISOString()
                  .slice(0, 10)}.json`}
                format="json"
                method="POST"
                downloadSize={summaryData?.total_case_count}
                filters={{
                  content: [
                    {
                      content: {
                        field: "files.file_id",
                        value: cart.map((file) => file.file_id),
                      },
                      op: "in",
                    },
                  ],
                  op: "and",
                }}
              />
              <Menu.Item
                component={DownloadButton}
                classNames={{ item: "font-normal" }}
                displayVariant="subtle"
                activeText="Processing"
                inactiveText="TSV"
                preventClickEvent
                showIcon={true}
                endpoint="biospecimen_tar"
                setActive={setBiospecimenTSVDownloadActive}
                active={biospecimenTSVDownloadActive}
                filename={`biospecimen.cart.${new Date()
                  .toISOString()
                  .slice(0, 10)}.tar.gz`}
                format="tsv"
                method="POST"
                downloadSize={summaryData?.total_case_count}
                filters={{
                  content: [
                    {
                      content: {
                        field: "files.file_id",
                        value: cart.map((file) => file.file_id),
                      },
                      op: "in",
                    },
                  ],
                  op: "and",
                }}
              />
            </Menu.Dropdown>
          </Menu>
          <Menu width="target">
            <Menu.Target>
              <Button
                data-testid="button-download-associated-data"
                classNames={{
                  root: `${buttonStyle} ${focusStyles}`,
                }}
                leftSection={
                  <DownloadIcon
                    aria-hidden="true"
                    size="1rem"
                    className="hidden xl:block"
                  />
                }
                rightSection={
                  <div className="border-l pl-1 -mr-2">
                    <ArrowDropDownIcon size="1.5em" aria-hidden="true" />
                  </div>
                }
              >
                Clinical
              </Button>
            </Menu.Target>
            <Menu.Dropdown data-testid="dropdown-menu-options">
              <Menu.Item
                component={DownloadButton}
                classNames={{ item: "font-normal" }}
                displayVariant="subtle"
                activeText="Processing"
                inactiveText="JSON"
                preventClickEvent
                showIcon={true}
                endpoint="clinical_tar"
                setActive={setClinicalJSONDownloadActive}
                active={clinicalJSONDownloadActive}
                filename={`clinical.cart.${new Date()
                  .toISOString()
                  .slice(0, 10)}.json`}
                format="json"
                method="POST"
                downloadSize={summaryData?.total_case_count}
                filters={{
                  content: [
                    {
                      content: {
                        field: "files.file_id",
                        value: cart.map((file) => file.file_id),
                      },
                      op: "in",
                    },
                  ],
                  op: "and",
                }}
              />
              <Menu.Item
                component={DownloadButton}
                classNames={{ item: "font-normal" }}
                displayVariant="subtle"
                activeText="Processing"
                inactiveText="TSV"
                preventClickEvent
                showIcon={true}
                endpoint="clinical_tar"
                setActive={setClinicalTSVDownloadActive}
                active={clinicalTSVDownloadActive}
                filename={`clinical.cart.${new Date()
                  .toISOString()
                  .slice(0, 10)}.tar.gz`}
                format="tsv"
                method="POST"
                downloadSize={summaryData?.total_case_count}
                filters={{
                  content: [
                    {
                      content: {
                        field: "files.file_id",
                        value: cart.map((file) => file.file_id),
                      },
                      op: "in",
                    },
                  ],
                  op: "and",
                }}
              />
            </Menu.Dropdown>
          </Menu>

          <DownloadButton
            displayVariant="header"
            activeText="Processing"
            inactiveText="Sample Sheet"
            preventClickEvent
            showIcon={true}
            endpoint="files"
            setActive={setSampleSheetDownloadActive}
            active={sampleSheetDownloadActive}
            filename={`gdc_sample_sheet.${new Date()
              .toISOString()
              .slice(0, 10)}.tsv`}
            format="tsv"
            method="POST"
            fields={[
              "file_id",
              "file_name",
              "data_category",
              "data_type",
              "cases.project.project_id",
              "cases.submitter_id",
              "cases.samples.submitter_id",
              "cases.samples.tissue_type",
              "cases.samples.tumor_descriptor",
              "cases.samples.specimen_type",
              "cases.samples.preservation_method",
            ]}
            filters={{
              content: [
                {
                  content: {
                    field: "files.file_id",
                    value: cart.map((file) => file.file_id),
                  },
                  op: "in",
                },
              ],
              op: "and",
            }}
            extraParams={{
              tsv_format: "sample-sheet",
            }}
          />
          <DownloadButton
            activeText="Processing"
            inactiveText="Metadata"
            showIcon={true}
            displayVariant="header"
            preventClickEvent
            endpoint="files"
            setActive={setMetadataDownloadActive}
            active={metadataDownloadActive}
            filename={`metadata.cart.${new Date()
              .toISOString()
              .slice(0, 10)}.json`}
            method="POST"
            filters={{
              content: [
                {
                  content: {
                    field: "files.file_id",
                    value: cart.map((file) => file.file_id),
                  },
                  op: "in",
                },
              ],
              op: "and",
            }}
            fields={[
              "state",
              "access",
              "md5sum",
              "data_format",
              "data_type",
              "data_category",
              "file_name",
              "file_size",
              "file_id",
              "platform",
              "experimental_strategy",
              "center.short_name",
              "annotations.annotation_id",
              "annotations.entity_id",
              "tags",
              "submitter_id",
              "archive.archive_id",
              "archive.submitter_id",
              "archive.revision",
              "associated_entities.entity_id",
              "associated_entities.entity_type",
              "associated_entities.case_id",
              "analysis.analysis_id",
              "analysis.workflow_type",
              "analysis.updated_datetime",
              "analysis.input_files.file_id",
              "analysis.metadata.read_groups.read_group_id",
              "analysis.metadata.read_groups.is_paired_end",
              "analysis.metadata.read_groups.read_length",
              "analysis.metadata.read_groups.library_name",
              "analysis.metadata.read_groups.sequencing_center",
              "analysis.metadata.read_groups.sequencing_date",
              "downstream_analyses.output_files.access",
              "downstream_analyses.output_files.file_id",
              "downstream_analyses.output_files.file_name",
              "downstream_analyses.output_files.data_category",
              "downstream_analyses.output_files.data_type",
              "downstream_analyses.output_files.data_format",
              "downstream_analyses.workflow_type",
              "downstream_analyses.output_files.file_size",
              "index_files.file_id",
            ]}
            extraParams={{
              expand: [
                "metadata_files",
                "annotations",
                "archive",
                "associated_entities",
                "center",
                "analysis",
                "analysis.input_files",
                "analysis.metadata",
                "analysis.metadata_files",
                "analysis.downstream_analyses",
                "analysis.downstream_analyses.output_files",
                "reference_genome",
                "index_file",
              ].join(","),
            }}
          />
          <Menu>
            <Menu.Target>
              <Button
                data-testid="button-remove-from-cart"
                leftSection={
                  <CartIcon
                    aria-hidden="true"
                    size="1rem"
                    className="hidden xl:block"
                  />
                }
                rightSection={
                  <div className="border-l pl-1 -mr-2">
                    <ArrowDropDownIcon size="1.5em" aria-hidden="true" />
                  </div>
                }
                classNames={{
                  root: `bg-nci-red-darker text-base-max hover:bg-removeButtonHover ${focusStyles}`,
                }}
              >
                Remove From Cart
              </Button>
            </Menu.Target>
            <Menu.Dropdown data-testid="dropdown-menu-options">
              <Menu.Item onClick={() => removeFromCart(cart, cart, dispatch)}>
                All Files ({cart.length})
              </Menu.Item>
              <Menu.Item
                onClick={() =>
                  removeFromCart(filesByCanAccess?.false || [], cart, dispatch)
                }
              >
                Unauthorized Files ({(filesByCanAccess?.false || []).length})
              </Menu.Item>
            </Menu.Dropdown>
          </Menu>
        </div>

        <h1 className="uppercase flex 2xl:ml-auto items-center truncate text-xl">
          Total of{" "}
          <FileIcon size={25} className="ml-2 mr-1" aria-hidden="true" />{" "}
          <b data-testid="text-file-count" className="mr-1">
            {summaryData?.total_doc_count?.toLocaleString() || "--"}
          </b>{" "}
          {summaryData?.total_doc_count === 1 ? "File" : "Files"}
          <PersonIcon size={25} className="ml-2 mr-1" aria-hidden="true" />{" "}
          <b data-testid="text-case-count" className="mr-1">
            {summaryData?.total_case_count?.toLocaleString() || "--"}
          </b>{" "}
          {summaryData?.total_case_count === 1 ? "Case" : "Cases"}{" "}
          <SaveIcon size={25} className="ml-2 mr-1" aria-hidden="true" />{" "}
          <span data-testid="text-size-count">
            {fileSize(summaryData?.total_file_size || 0)}
          </span>{" "}
        </h1>
      </div>
    </>
  );
};

export default CartHeader;
