import Link from "next/link";
import {
  AnnotationDefaults,
  CartFile,
  Union,
  useCoreDispatch,
} from "@gff/core";
import { ColumnDef, createColumnHelper } from "@tanstack/react-table";
import { Dispatch, SetStateAction, useMemo } from "react";
import { Button, Checkbox, Menu, Tooltip } from "@mantine/core";
import { allFilesInCart } from "@/utils/index";
import { addToCart, removeFromCart } from "@/features/cart/updateCart";
import { ImageSlideCount } from "@/components/ImageSlideCount";
import OverflowTooltippedLabel from "@/components/OverflowTooltippedLabel";
import { PopupIconButton } from "@/components/PopupIconButton/PopupIconButton";
import { entityMetadataType } from "@/utils/contexts";
import { ArraySeparatedSpan } from "@/components/ArraySeparatedSpan/ArraySeparatedSpan";
import {
  CartIcon,
  DropdownIcon,
  TrashIcon,
  AddToQueueIcon,
} from "@/utils/icons";

export type casesTableDataType = {
  case_uuid: string;
  case_id: string;
  project: string;
  program: string;
  primary_site: string;
  disease_type: string;
  primary_diagnosis: string;
  age_at_diagnosis: string;
  vital_status: string;
  days_to_death: string;
  gender: string;
  race: string;
  ethnicity: string;
  slide_count: number;
  files_count: number;
  files: {
    access: "open" | "controlled";
    acl: string[];
    file_id: string;
    file_size: number;
    state: string;
    file_name: string;
    data_type: string;
  }[];
  experimental_strategies: (string | number)[] | string;
  annotations: AnnotationDefaults[];
};

const casesDataColumnHelper = createColumnHelper<casesTableDataType>();

export const useGenerateCasesTableColumns = ({
  currentCart,
  setEntityMetadata,
  currentPage,
  totalPages,
}: {
  currentCart: CartFile[];
  setEntityMetadata: Dispatch<SetStateAction<entityMetadataType>>;
  currentPage: number;
  totalPages: number;
}): ColumnDef<casesTableDataType>[] => {
  const dispatch = useCoreDispatch();
  const CasesTableDefaultColumns = useMemo<ColumnDef<casesTableDataType>[]>(
    () => [
      casesDataColumnHelper.display({
        id: "select",
        header: ({ table }) => (
          <Checkbox
            data-testid="checkbox-select-all-cases-table"
            size="xs"
            classNames={{
              input: "checked:bg-accent checked:border-accent",
              label: "sr-only",
            }}
            label={`Select all the case rows on page ${currentPage} of ${totalPages}`}
            {...{
              checked: table.getIsAllRowsSelected(),
              onChange: table.getToggleAllRowsSelectedHandler(),
            }}
          />
        ),
        cell: ({ row }) => (
          <Checkbox
            data-testid="checkbox-select-cases-table"
            size="xs"
            classNames={{
              input: "checked:bg-accent checked:border-accent",
            }}
            aria-label={row?.original?.project}
            {...{
              checked: row.getIsSelected(),
              onChange: row.getToggleSelectedHandler(),
            }}
          />
        ),
      }),
      casesDataColumnHelper.display({
        id: "cart",
        header: "Cart",
        cell: ({ row }) => {
          const isAllFilesInCart = row.original.files
            ? allFilesInCart(currentCart, row.original.files)
            : false;
          const numberOfFilesToRemove = row.original.files
            ?.map((file) =>
              currentCart.filter(
                (data_file) => data_file.file_id === file.file_id,
              ),
            )
            .filter((item) => item.length > 0).length;
          const isDisabled = row.original.files_count === 0;
          const isPlural = row.original.files_count > 1;
          return (
            <Menu position="bottom-start" zIndex={300}>
              <Menu.Target>
                <Tooltip label="No files to add to Cart" disabled={!isDisabled}>
                  <Button
                    data-testid="button-add-remove-cases-table"
                    aria-label={`${
                      isAllFilesInCart ? "remove" : "add"
                    } all files ${isAllFilesInCart ? "from" : "to"} the cart`}
                    leftSection={
                      <div className="mr-2">
                        <CartIcon
                          className={
                            isAllFilesInCart && "text-primary-contrast-darkest"
                          }
                          aria-hidden="true"
                        />
                      </div>
                    }
                    rightSection={
                      <div className="border-l">
                        <DropdownIcon
                          className={
                            isAllFilesInCart && "text-primary-contrast-darkest"
                          }
                          size={18}
                          aria-hidden="true"
                        />
                      </div>
                    }
                    variant="outline"
                    classNames={{
                      root: "w-12 pr-0 bg-base-max text-primary disabled:border disabled:bg-base-lightest disabled:opacity-50 disabled:border-primary",
                      section: "m-0",
                    }}
                    size="compact-xs"
                    className={`${isAllFilesInCart && "bg-primary-darkest"}`}
                    disabled={isDisabled}
                  />
                </Tooltip>
              </Menu.Target>
              <Menu.Dropdown>
                {numberOfFilesToRemove < row.original.files_count && (
                  <Menu.Item
                    data-testid="button-add-files-to-cart-cases-table"
                    leftSection={<AddToQueueIcon />}
                    onClick={() => {
                      addToCart(row.original.files, currentCart, dispatch);
                    }}
                  >
                    Add {row.original.files_count - numberOfFilesToRemove} Case{" "}
                    {isPlural ? "files" : "file"} to the Cart
                  </Menu.Item>
                )}

                {numberOfFilesToRemove > 0 && (
                  <Menu.Item
                    data-testid="button-remove-files-from-cart-cases-table"
                    leftSection={<TrashIcon />}
                    onClick={() => {
                      removeFromCart(row.original.files, currentCart, dispatch);
                    }}
                  >
                    Remove{" "}
                    {numberOfFilesToRemove === 0
                      ? row.original.files_count
                      : numberOfFilesToRemove}{" "}
                    Case {isPlural ? "files" : "file"} from the Cart
                  </Menu.Item>
                )}
              </Menu.Dropdown>
            </Menu>
          );
        },
      }),
      casesDataColumnHelper.display({
        id: "slides",
        header: "Slides",
        cell: ({ row }) => (
          // This needs both passHref and legacyBehavior: https://nextjs.org/docs/pages/api-reference/components/link#if-the-child-is-a-functional-component
          <Link
            href={{
              pathname: "/image-viewer/MultipleImageViewerPage",
              query: { caseId: row.original.case_uuid },
            }}
            passHref
            legacyBehavior
          >
            <ImageSlideCount slideCount={row.original.slide_count} />
          </Link>
        ),
      }),
      casesDataColumnHelper.accessor("case_id", {
        id: "case_id",
        header: "Case ID",
        cell: ({ row }) => (
          <OverflowTooltippedLabel label={row.original.case_id}>
            <PopupIconButton
              handleClick={() =>
                setEntityMetadata({
                  entity_type: "case",
                  entity_id: row.original.case_uuid,
                })
              }
              label={row.original.case_id}
            />
          </OverflowTooltippedLabel>
        ),
      }),
      casesDataColumnHelper.accessor("case_uuid", {
        id: "case_uuid",
        header: "Case UUID",
      }),
      casesDataColumnHelper.accessor("project", {
        id: "project",
        header: "Project",
        cell: ({ row }) => (
          <OverflowTooltippedLabel label={row.original.project}>
            <PopupIconButton
              handleClick={() =>
                setEntityMetadata({
                  entity_type: "project",
                  entity_id: row.original.project,
                })
              }
              label={row.original.project}
            />
          </OverflowTooltippedLabel>
        ),
      }),
      casesDataColumnHelper.accessor("program", {
        id: "program",
        header: "Program",
      }),
      casesDataColumnHelper.accessor("primary_site", {
        id: "primary_site",
        header: "Primary Site",
      }),
      casesDataColumnHelper.accessor("disease_type", {
        id: "disease_type",
        header: "Disease Type",
      }),
      casesDataColumnHelper.accessor("primary_diagnosis", {
        id: "primary_diagnosis",
        header: "Primary Diagnosis",
        enableSorting: false,
      }),
      casesDataColumnHelper.accessor("age_at_diagnosis", {
        id: "age_at_diagnosis",
        header: "Age At Diagnosis",
        enableSorting: false,
      }),
      casesDataColumnHelper.accessor("vital_status", {
        id: "vital_status",
        header: "Vital Status",
      }),
      casesDataColumnHelper.accessor("days_to_death", {
        id: "days_to_death",
        header: "Days To Death",
      }),
      casesDataColumnHelper.accessor("gender", {
        id: "gender",
        header: "Gender",
      }),
      casesDataColumnHelper.accessor("race", {
        id: "race",
        header: "Race",
      }),
      casesDataColumnHelper.accessor("ethnicity", {
        id: "ethnicity",
        header: "Ethnicity",
      }),
      casesDataColumnHelper.accessor("files_count", {
        id: "files",
        header: "Files",
        cell: ({ getValue }) => getValue().toLocaleString(),
      }),
      casesDataColumnHelper.accessor("experimental_strategies", {
        id: "experimental_strategy",
        header: "Experimental Strategy",
        cell: ({ getValue }) => (
          <ArraySeparatedSpan data={getValue() as string[]} />
        ),
        enableSorting: false,
      }),
      casesDataColumnHelper.display({
        id: "annotations",
        header: "Annotations",
        cell: ({ row }) => row.original.annotations.length,
      }),
    ],
    [currentCart, dispatch, setEntityMetadata, currentPage, totalPages],
  );

  return CasesTableDefaultColumns;
};

export const buildCasesTableSearchFilters = (
  term?: string,
): Union | undefined => {
  if (term !== undefined && term.length > 0) {
    return {
      operator: "or",
      operands: [
        {
          operator: "includes",
          field: "cases.case_id", // case insensitive
          operands: [`*${term}*`],
        },
        {
          operator: "includes",
          field: "cases.submitter_id", // case sensitive
          operands: [`*${term}*`],
        },
      ],
    };
  }
  return undefined;
};

export const MAX_CASE_IDS = 100000;
