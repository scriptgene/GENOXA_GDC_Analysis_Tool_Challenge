import { createAsyncThunk, createSlice } from "@reduxjs/toolkit";
import {
  CoreDataSelectorResponse,
  createUseFiltersCoreDataHook,
  DataStatus,
} from "../../dataAccess";
import { castDraft } from "immer";
import { CoreDispatch } from "../../store";
import { CoreState } from "../../reducers";
import { fetchSmsAggregations } from "./smsAggregationsApi";
import { GraphQLApiResponse, graphqlAPI } from "../gdcapi/gdcgraphql";
import { GenomicTableProps } from "./types";
import {
  buildCohortGqlOperator,
  extractFiltersWithPrefixFromFilterSet,
  filterSetToOperation,
  selectCurrentCohortFilters,
} from "../cohort";
import {
  convertFilterToGqlFilter,
  Union,
  UnionOrIntersection,
} from "../gdcapi/filters";

const GenesTableGraphQLQuery = `
          query GenesTable(
            $caseFilters: FiltersArgument
            $genesTable_filters: FiltersArgument
            $genesTable_size: Int
            $genesTable_offset: Int
            $score: String
            $ssmCase: FiltersArgument
            $geneCaseFilter: FiltersArgument
            $ssmTested: FiltersArgument
            $cnvTested: FiltersArgument
            $cnvGainFilters: FiltersArgument
            $cnvLossFilters: FiltersArgument
            $cnvAmplificationFilters: FiltersArgument
            $cnvHomozygousDeletionFilters: FiltersArgument
            $sort: [Sort]
          ) {
            genesTableViewer: viewer {
              explore {
                cases {
                  hits(first: 0, case_filters: $ssmTested) {
                    total
                  }
                }
                filteredCases: cases {
                  hits(first: 0, case_filters: $geneCaseFilter) {
                    total
                  }
                }
                cnvCases: cases {
                  hits(first: 0, case_filters: $cnvTested) {
                    total
                  }
                }
                genes {
                  hits(
                    first: $genesTable_size
                    offset: $genesTable_offset
                    filters: $genesTable_filters
                    case_filters: $caseFilters
                    score: $score
                    sort: $sort
                  ) {
                    total
                    edges {
                      node {
                        id
                        numCases: score
                        symbol
                        name
                        cytoband
                        biotype
                        gene_id
                        is_cancer_gene_census
                        ssm_case: case {
                          hits(first: 0, filters: $ssmCase) {
                            total
                          }
                        }
                        cnv_case: case {
                          hits(first: 0, case_filters: $caseFilters, filters: $cnvTested) {
                            total
                          }
                        }
                        case_cnv_gain: case {
                          hits(first: 0, case_filters: $caseFilters, filters: $cnvGainFilters) {
                            total
                          }
                        }
                        case_cnv_loss: case {
                          hits(first: 0, case_filters: $caseFilters, filters: $cnvLossFilters) {
                            total
                          }
                        }
                        case_cnv_amplification: case {
                          hits(first: 0, case_filters: $caseFilters, filters: $cnvAmplificationFilters) {
                            total
                          }
                        }
                        case_cnv_homozygous_deletion: case {
                          hits(first: 0, case_filters: $caseFilters, filters: $cnvHomozygousDeletionFilters) {
                            total
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
`;

export interface GeneRowInfo {
  readonly biotype: string;
  readonly case_cnv_amplification: number;
  readonly case_cnv_gain: number;
  readonly case_cnv_loss: number;
  readonly case_cnv_homozygous_deletion: number;
  readonly cnv_case: number;
  readonly cytoband: Array<string>;
  readonly gene_id: string;
  readonly id: string;
  readonly is_cancer_gene_census: boolean;
  readonly name: string;
  readonly numCases: number;
  readonly ssm_case: number;
  readonly symbol: string;
}

export interface GDCGenesTable {
  readonly cases: number;
  readonly cnvCases: number;
  readonly filteredCases: number;
  readonly genes: ReadonlyArray<GeneRowInfo>;
  readonly genes_total: number;
  readonly mutationCounts?: Record<string, string>;
}

export type CnvChange =
  | "Amplification"
  | "Gain"
  | "Loss"
  | "Homozygous Deletion";

export const buildGeneTableSearchFilters = (
  term?: string,
): Union | undefined => {
  if (term !== undefined) {
    return {
      operator: "or",
      operands: [
        {
          operator: "includes",
          field: "genes.cytoband",
          operands: [`*${term}*`],
        },
        {
          operator: "includes",
          field: "genes.gene_id",
          operands: [`*${term}*`],
        },
        {
          operator: "includes",
          field: "genes.symbol",
          operands: [`*${term}*`],
        },
        { operator: "includes", field: "genes.name", operands: [`*${term}*`] },
      ],
    };
  }
  return undefined;
};

export const fetchGenesTable = createAsyncThunk<
  GraphQLApiResponse,
  GenomicTableProps,
  { dispatch: CoreDispatch; state: CoreState }
>(
  "genes/genesTable",
  async (
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    {
      pageSize,
      offset,
      genesTableFilters,
      genomicFilters,
      cohortFilters,
    }: GenomicTableProps,
  ): Promise<GraphQLApiResponse> => {
    const caseFilters = buildCohortGqlOperator(cohortFilters);
    const cohortFiltersContent = caseFilters?.content
      ? Object(caseFilters?.content)
      : [];
    const genesTable_filters = buildCohortGqlOperator(genesTableFilters);

    const baseFilters = filterSetToOperation(genomicFilters) as
      | UnionOrIntersection
      | undefined;

    const rawFilterContents =
      baseFilters && convertFilterToGqlFilter(baseFilters)?.content;
    const filterContents = rawFilterContents ? Object(rawFilterContents) : [];

    /**
     * Only apply "genes." filters to the genes table's CNV gain and loss filters.
     */
    const onlyGenesFilters = buildCohortGqlOperator(
      extractFiltersWithPrefixFromFilterSet(genomicFilters, "genes."),
    );

    const graphQlFilters = {
      caseFilters: caseFilters ? caseFilters : {},
      genesTable_filters: genesTable_filters ? genesTable_filters : {},
      genesTable_size: pageSize,
      genesTable_offset: offset,
      score: "case.project.project_id",
      ssmCase: {
        op: "and",
        content: [
          {
            op: "in",
            content: {
              field: "cases.available_variation_data",
              value: ["ssm"],
            },
          },
          {
            op: "NOT",
            content: {
              field: "genes.case.ssm.observation.observation_id",
              value: "MISSING",
            },
          },
        ],
      },
      geneCaseFilter: {
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["ssm"],
              },
              op: "in",
            },
          ],
          ...cohortFiltersContent,
        ],
        op: "and",
      },
      ssmTested: {
        content: [
          {
            content: {
              field: "cases.available_variation_data",
              value: ["ssm"],
            },
            op: "in",
          },
        ],
        op: "and",
      },
      cnvTested: {
        op: "and",
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["cnv"],
              },
              op: "in",
            },
          ],
          ...cohortFiltersContent,
        ],
      },
      cnvGainFilters: {
        op: "and",
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["cnv"],
              },
              op: "in",
            },
            {
              content: {
                field: "cnvs.cnv_change_5_category",
                value: ["Gain"],
              },
              op: "in",
            },
          ],
          ...(onlyGenesFilters?.content
            ? Object(onlyGenesFilters?.content)
            : []),
        ],
      },
      cnvLossFilters: {
        op: "and",
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["cnv"],
              },
              op: "in",
            },
            {
              content: {
                field: "cnvs.cnv_change_5_category",
                value: ["Loss"],
              },
              op: "in",
            },
          ],
          ...(onlyGenesFilters?.content
            ? Object(onlyGenesFilters?.content)
            : []),
        ],
      },
      cnvAmplificationFilters: {
        op: "and",
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["cnv"],
              },
              op: "in",
            },
            {
              content: {
                field: "cnvs.cnv_change_5_category",
                value: ["Amplification"],
              },
              op: "in",
            },
          ],
          ...(onlyGenesFilters?.content
            ? Object(onlyGenesFilters?.content)
            : []),
        ],
      },
      cnvHomozygousDeletionFilters: {
        op: "and",
        content: [
          ...[
            {
              content: {
                field: "cases.available_variation_data",
                value: ["cnv"],
              },
              op: "in",
            },
            {
              content: {
                field: "cnvs.cnv_change_5_category",
                value: ["Homozygous Deletion"],
              },
              op: "in",
            },
          ],
          ...(onlyGenesFilters?.content
            ? Object(onlyGenesFilters?.content)
            : []),
        ],
      },
    };

    // get the TableData

    const results: GraphQLApiResponse<any> = await graphqlAPI(
      GenesTableGraphQLQuery,
      graphQlFilters,
    );
    // if we have valid data from the table, need to query the ssms counts
    if (!results.errors) {
      // extract the gene ids and user it for the call to
      const geneIds =
        results.data.genesTableViewer.explore.genes.hits.edges.map(
          ({ node }: Record<string, any>) => node.gene_id,
        );
      const counts = await fetchSmsAggregations({
        ids: geneIds,
        filters: filterContents,
        caseFilters: caseFilters,
      });
      if (!counts.errors) {
        const countsData =
          counts.data.ssmsAggregationsViewer.explore.ssms.aggregations
            .consequence__transcript__gene__gene_id;
        results.data.genesTableViewer["mutationCounts"] =
          countsData.buckets.reduce(
            (
              counts: Record<string, number>,
              apiBucket: Record<string, any>,
            ) => {
              counts[apiBucket.key] = apiBucket.doc_count.toLocaleString();
              return counts;
            },
            {} as Record<string, number>,
          );
      }
    }

    return results;
  },
);

export interface GenesTableState {
  readonly genes: GDCGenesTable;
  readonly status: DataStatus;
  readonly error?: string;
  readonly requestId?: string;
}

const initialState: GenesTableState = {
  genes: {
    cases: 0,
    filteredCases: 0,
    cnvCases: 0,
    genes: [],
    genes_total: 0,
  },
  status: "uninitialized",
};

const slice = createSlice({
  name: "genomic/genesTable",
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(fetchGenesTable.fulfilled, (state, action) => {
        if (state.requestId != action.meta.requestId) return state;
        const response = action.payload;
        if (response.errors) {
          state = castDraft(initialState);
          state.status = "rejected";
          state.error = response.errors.filters;
        }

        const data = action.payload.data.genesTableViewer.explore;
        state.genes.cases = data.cases.hits.total;
        state.genes.cnvCases = data.cnvCases.hits.total;
        state.genes.filteredCases = data.filteredCases.hits.total;
        state.genes.genes_total = data.genes.hits.total;
        state.genes.mutationCounts =
          action.payload.data.genesTableViewer.mutationCounts;
        state.genes.genes = data.genes.hits.edges.map(
          ({ node }: Record<string, any>): GeneRowInfo => {
            const {
              biotype,
              cytoband,
              gene_id,
              id,
              is_cancer_gene_census,
              name,
              numCases,
              symbol,
            } = node;
            return {
              biotype,
              cytoband,
              gene_id,
              id,
              is_cancer_gene_census,
              name,
              numCases,
              symbol,
              cnv_case: node.cnv_case.hits.total,
              case_cnv_amplification: node.case_cnv_amplification.hits.total,
              case_cnv_loss: node.case_cnv_loss.hits.total,
              case_cnv_gain: node.case_cnv_gain.hits.total,
              case_cnv_homozygous_deletion:
                node.case_cnv_homozygous_deletion.hits.total,
              ssm_case: node.ssm_case.hits.total,
            };
          },
        );
        state.status = "fulfilled";
        state.error = undefined;
        return state;
      })
      .addCase(fetchGenesTable.pending, (state, action) => {
        state.status = "pending";
        state.requestId = action.meta.requestId;
        return state;
      })
      .addCase(fetchGenesTable.rejected, (state, action) => {
        if (state.requestId != action.meta.requestId) return state;
        state.status = "rejected";
        if (action.error) {
          state.error = action.error.message;
        }
        return state;
      });
  },
});

export const genesTableReducer = slice.reducer;

export const selectGenesTableState = (state: CoreState): GDCGenesTable =>
  state.genomic.genesTable.genes;

export const selectGenesTableData = (
  state: CoreState,
): CoreDataSelectorResponse<GenesTableState> => {
  return {
    data: state.genomic.genesTable,
    status: state.genomic.genesTable.status,
    error: state.genomic.genesTable.error,
  };
};

export const useGenesTable = createUseFiltersCoreDataHook(
  fetchGenesTable,
  selectGenesTableData,
  selectCurrentCohortFilters,
);
