import { Middleware, Reducer } from "@reduxjs/toolkit";
import { castDraft } from "immer";
import { endpointSlice } from "../gdcapi";
import { GdcApiRequest, GdcApiResponse } from "../gdcapi/types";
import { FileDefaults } from "../gdcapi/types";

const accessTypes = ["open", "controlled"] as const;

export type AccessType = (typeof accessTypes)[number];

const isAccessType = (x: unknown): x is AccessType => {
  return accessTypes.some((t) => t === x);
};

const asAccessType = (x: unknown): AccessType => {
  if (isAccessType(x)) {
    return x;
  } else {
    throw new Error(`${x} is not a valid access type`);
  }
};

// TODO use CartFile instead and combine anything that submits to cart
export interface GdcCartFile {
  readonly file_name: string;
  readonly data_category: string;
  readonly data_type: string;
  readonly data_format: string;
  readonly state: string;
  readonly file_size: number;
  readonly file_id: string;
  readonly access: AccessType;
  readonly acl: ReadonlyArray<string>;
  readonly project_id?: string;
  readonly createdDatetime: string;
  readonly updatedDatetime: string;
  readonly submitterId: string;
  readonly md5sum: string;
}
export type FileAnnotationsType = {
  readonly annotation_id: string;
  readonly case_id?: string;
  readonly case_submitter_id?: string;
  readonly category: string;
  readonly classification: string;
  readonly created_datetime: string;
  readonly entity_id: string;
  readonly entity_submitter_id: string;
  readonly entity_type: string;
  readonly notes: string;
  readonly state: string;
  readonly status: string;
  readonly updated_datetime: string;
};
export type FileCaseType = ReadonlyArray<{
  readonly case_id: string;
  readonly submitter_id: string;
  readonly annotations?: ReadonlyArray<string>;
  readonly project?: {
    readonly dbgap_accession_number?: string;
    readonly disease_type: string;
    readonly name: string;
    readonly primary_site: string;
    readonly project_id: string;
    readonly releasable: boolean;
    readonly released: boolean;
    readonly state: string;
  };
  readonly samples?: ReadonlyArray<{
    readonly sample_id: string;
    readonly submitter_id: string;
    readonly tissue_type: string;
    readonly tumor_descriptor: string;
    readonly portions?: ReadonlyArray<{
      readonly submitter_id: string;
      readonly analytes?: ReadonlyArray<{
        readonly analyte_id: string;
        readonly analyte_type: string;
        readonly submitter_id: string;
        readonly aliquots?: ReadonlyArray<{
          readonly aliquot_id: string;
          readonly submitter_id: string;
        }>;
      }>;
      readonly slides?: ReadonlyArray<{
        readonly created_datetime: string | null;
        readonly number_proliferating_cells: number | null;
        readonly percent_eosinophil_infiltration: number | null;
        readonly percent_granulocyte_infiltration: number | null;
        readonly percent_inflam_infiltration: number | null;
        readonly percent_lymphocyte_infiltration: number | null;
        readonly percent_monocyte_infiltration: number | null;
        readonly percent_neutrophil_infiltration: number | null;
        readonly percent_necrosis: number | null;
        readonly percent_normal_cells: number | null;
        readonly percent_stromal_cells: number | null;
        readonly percent_tumor_cells: number | null;
        readonly percent_tumor_nuclei: number | null;
        readonly section_location: string | null;
        readonly slide_id: string | null;
        readonly state: string | null;
        readonly submitter_id: string | null;
        readonly updated_datetime: string | null;
      }>;
    }>;
  }>;
}>;
export interface GdcFile {
  readonly id?: string;
  readonly submitterId: string;
  readonly access: AccessType;
  readonly acl: ReadonlyArray<string>;
  readonly createdDatetime: string;
  readonly updatedDatetime: string;
  readonly data_category: string;
  readonly data_format: string;
  readonly dataRelease?: string;
  readonly data_type: string;
  readonly file_id: string;
  readonly file_name: string;
  readonly file_size: number;
  readonly md5sum: string;
  readonly platform?: string;
  readonly state: string;
  readonly fileType?: string;
  readonly version?: string;
  readonly experimental_strategy?: string;
  readonly project_id?: string;
  readonly annotations?: ReadonlyArray<FileAnnotationsType>;
  readonly cases?: FileCaseType;
  readonly total_reads?: string;
  readonly average_base_quality?: number;
  readonly average_insert_size?: number;
  readonly average_read_length?: number;
  readonly mean_coverage?: number;
  readonly pairs_on_diff_chr?: string;
  readonly contamination?: number;
  readonly contamination_error?: number;
  readonly proportion_reads_mapped?: number;
  readonly proportion_reads_duplicated?: number;
  readonly proportion_base_mismatch?: number;
  readonly proportion_targets_no_coverage?: number;
  readonly proportion_coverage_10x?: number;
  readonly proportion_coverage_30x?: number;
  readonly msi_score?: number;
  readonly msi_status?: number;
  readonly associated_entities?: ReadonlyArray<{
    readonly entity_submitter_id: string;
    readonly entity_type: string;
    readonly case_id: string;
    readonly entity_id: string;
  }>;
  readonly analysis?: {
    readonly workflow_type: string;
    readonly input_files?: GdcCartFile[];
    readonly metadata?: {
      readonly read_groups: Array<{
        readonly read_group_id: string;
        readonly is_paired_end: boolean;
        readonly read_length: number;
        readonly library_name: string;
        readonly sequencing_center: string;
        readonly sequencing_date: string;
      }>;
    };
  };
  readonly downstream_analyses?: ReadonlyArray<{
    readonly workflow_type: string;
    readonly output_files?: GdcCartFile[];
  }>;
  readonly index_files?: ReadonlyArray<{
    readonly submitterId: string;
    readonly createdDatetime: string;
    readonly updatedDatetime: string;
    readonly data_category: string;
    readonly data_format: string;
    readonly data_type: string;
    readonly file_id: string;
    readonly file_name: string;
    readonly file_size: number;
    readonly md5sum: string;
    readonly state: string;
  }>;
}

export const mapFileData = (files: ReadonlyArray<FileDefaults>): GdcFile[] => {
  return files.map((hit) => ({
    id: hit.id,
    submitterId: hit.submitter_id,
    access: asAccessType(hit.access),
    acl: [...hit.acl],
    createdDatetime: hit.created_datetime,
    updatedDatetime: hit.updated_datetime,
    data_category: hit.data_category,
    data_format: hit.data_format,
    dataRelease: hit.data_release,
    data_type: hit.data_type,
    file_id: hit.file_id,
    file_name: hit.file_name,
    file_size: hit.file_size,
    md5sum: hit.md5sum,
    platform: hit.platform,
    state: hit.state,
    fileType: hit.type,
    version: hit.version,
    experimental_strategy: hit.experimental_strategy,
    project_id: hit.cases?.[0]?.project?.project_id,
    annotations: hit.annotations?.map((annotation) => annotation),
    total_reads: hit.total_reads?.toLocaleString(),
    average_insert_size: hit?.average_insert_size,
    average_read_length: hit?.average_read_length,
    average_base_quality: hit?.average_base_quality,
    mean_coverage: hit?.mean_coverage,
    pairs_on_diff_chr: hit?.pairs_on_diff_chr?.toLocaleString(),
    contamination: hit?.contamination,
    contamination_error: hit?.contamination_error,
    proportion_reads_mapped: hit?.proportion_reads_mapped,
    proportion_reads_duplicated: hit?.proportion_reads_duplicated,
    proportion_base_mismatch: hit?.proportion_base_mismatch,
    proportion_targets_no_coverage: hit?.proportion_targets_no_coverage,
    proportion_coverage_10x: hit?.proportion_coverage_10x,
    proportion_coverage_30x: hit?.proportion_coverage_30x,
    msi_score: hit?.msi_score,
    msi_status: hit?.msi_status,
    cases: hit.cases?.map((caseObj) => {
      return {
        case_id: caseObj.case_id,
        submitter_id: caseObj.submitter_id,
        annotations: caseObj.annotations?.map(
          (annotation) => annotation.annotation_id,
        ),
        project: caseObj?.project
          ? {
              dbgap_accession_number: caseObj.project.dbgap_accession_number,
              disease_type: caseObj.project.disease_type,
              name: caseObj.project.name,
              primary_site: caseObj.project.primary_site,
              project_id: caseObj.project.project_id,
              releasable: caseObj.project.releasable,
              released: caseObj.project.released,
              state: caseObj.project.state,
            }
          : undefined,
        samples: caseObj.samples?.map((sample) => {
          return {
            sample_id: sample.sample_id,
            submitter_id: sample.submitter_id,
            tissue_type: sample.tissue_type,
            tumor_descriptor: sample.tumor_descriptor,
            portions: sample.portions?.map((portion) => {
              return {
                submitter_id: portion.submitter_id,
                analytes: portion.analytes?.map((analyte) => {
                  return {
                    analyte_id: analyte.analyte_id,
                    analyte_type: analyte.analyte_type,
                    submitter_id: analyte.submitter_id,
                    aliquots: analyte.aliquots?.map((aliquot) => {
                      return {
                        aliquot_id: aliquot.aliquot_id,
                        submitter_id: aliquot.submitter_id,
                      };
                    }),
                  };
                }),
                slides: portion.slides?.map((slide) => {
                  return {
                    number_proliferating_cells:
                      slide.number_proliferating_cells,
                    percent_eosinophil_infiltration:
                      slide.percent_eosinophil_infiltration,
                    percent_granulocyte_infiltration:
                      slide.percent_granulocyte_infiltration,
                    percent_inflam_infiltration:
                      slide.percent_inflam_infiltration,
                    percent_lymphocyte_infiltration:
                      slide.percent_lymphocyte_infiltration,
                    percent_monocyte_infiltration:
                      slide.percent_monocyte_infiltration,
                    percent_necrosis: slide.percent_necrosis,
                    percent_neutrophil_infiltration:
                      slide.percent_neutrophil_infiltration,
                    percent_normal_cells: slide.percent_normal_cells,
                    percent_stromal_cells: slide.percent_stromal_cells,
                    percent_tumor_cells: slide.percent_tumor_cells,
                    percent_tumor_nuclei: slide.percent_tumor_nuclei,
                    section_location: slide.section_location,
                    slide_id: slide.slide_id,
                    state: slide.state,
                    submitter_id: slide.submitter_id,
                    updated_datetime: slide.updated_datetime,
                    created_datetime: slide.created_datetime,
                  };
                }),
              };
            }),
          };
        }),
      };
    }),
    associated_entities: hit.associated_entities?.map((associated_entity) => ({
      entity_submitter_id: associated_entity.entity_submitter_id,
      entity_type: associated_entity.entity_type,
      case_id: associated_entity.case_id,
      entity_id: associated_entity.entity_id,
    })),
    analysis: hit.analysis
      ? {
          workflow_type: hit.analysis.workflow_type,
          input_files: hit.analysis.input_files?.map((file) => {
            return {
              file_name: file.file_name,
              data_category: file.data_category,
              data_type: file.data_type,
              data_format: file.data_format,
              file_size: file.file_size,
              file_id: file.file_id,
              acl: hit.acl,
              access: asAccessType(file.access),
              project_id: hit.cases?.[0].project?.project_id,
              state: file.state,
              submitterId: file.submitter_id,
              createdDatetime: file.created_datetime,
              updatedDatetime: file.updated_datetime,
              md5sum: file.md5sum,
            };
          }),
          metadata: hit.analysis.metadata
            ? {
                read_groups: hit.analysis.metadata.read_groups.map(
                  (read_group) => ({
                    read_group_id: read_group.read_group_id,
                    is_paired_end: read_group.is_paired_end,
                    read_length: read_group.read_length,
                    library_name: read_group.library_name,
                    sequencing_center: read_group.sequencing_center,
                    sequencing_date: read_group.sequencing_date,
                  }),
                ),
              }
            : undefined,
        }
      : undefined,
    downstream_analyses: hit.downstream_analyses?.map((analysis) => {
      return {
        workflow_type: analysis.workflow_type,
        output_files: analysis.output_files?.map((file) => {
          return {
            file_name: file.file_name,
            data_category: file.data_category,
            data_type: file.data_type,
            data_format: file.data_format,
            file_size: file.file_size,
            file_id: file.file_id,
            acl: hit.acl,
            access: asAccessType(file.access),
            project_id: hit.cases?.[0].project?.project_id,
            state: file.state,
            submitterId: file.submitter_id,
            createdDatetime: file.created_datetime,
            updatedDatetime: file.updated_datetime,
            md5sum: file.md5sum,
          };
        }),
      };
    }),
    index_files: hit.index_files?.map((idx) => ({
      submitterId: idx.submitter_id,
      createdDatetime: idx.created_datetime,
      updatedDatetime: idx.updated_datetime,
      data_category: idx.data_category,
      data_format: idx.data_format,
      data_type: idx.data_type,
      file_id: idx.file_id,
      file_name: idx.file_name,
      file_size: idx.file_size,
      md5sum: idx.md5sum,
      state: idx.state,
    })),
  }));
};

export const filesApiSlice = endpointSlice.injectEndpoints({
  endpoints: (builder) => ({
    getFiles: builder.query({
      query: (request: GdcApiRequest) => ({
        request,
        endpoint: "files",
        fetchAll: false,
      }),
      transformResponse: (response: GdcApiResponse<FileDefaults>) => {
        if (response.warnings && Object.keys(response.warnings).length > 0)
          return {
            files: [],
            pagination: undefined,
          };

        return {
          files: castDraft(mapFileData(response?.data?.hits ?? [])),
          pagination: response.data.pagination,
        };
      },
    }),
  }),
});
export const { useGetFilesQuery } = filesApiSlice;

export const filesApiSliceMiddleware = filesApiSlice.middleware as Middleware;
export const filesApiSliceReducerPath: string = filesApiSlice.reducerPath;
export const filesApiReducer: Reducer = filesApiSlice.reducer as Reducer;
