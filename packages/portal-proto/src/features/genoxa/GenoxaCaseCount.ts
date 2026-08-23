import { buildCohortGqlOperator, FilterSet, graphqlAPISlice } from "@gff/core";

/**
 * Single source of truth for the GENOXA prefilter.
 * Import this in the GENOXA app query too, so the card count
 * and the analysis always describe the same set of cases.
 */
export const GENOXA_FILE_FILTER = {
  op: "and",
  content: [
    {
        "op": "in",
        "content": {
          "field": "cases.samples.tissue_type",
          "value": ["tumor","normal"]
        }
    },
    
    {
      op: "in",
      content: {
        field: "files.data_category",
        value: ["transcriptome profiling"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.data_type",
        value: ["Gene Expression Quantification"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.experimental_strategy",
        value: ["RNA-Seq"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.analysis.workflow_type",
        value: ["STAR - Counts"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.data_format",
        value: ["tsv"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.platform",
        value: ["illumina"],
      },
    },
    {
      op: "in",
      content: {
        field: "files.access",
        value: ["open"],
      },
    },
  ],
};

const graphQLQuery = `
  query genoxaCaseCountQuery($cohortFilters: FiltersArgument,
  $genoxaCaseFilter: FiltersArgument) {
  viewer {
    explore {
      genoxaCaseCount : cases {
        hits(case_filters: $cohortFilters, filters: $genoxaCaseFilter, first: 0) {
          total
        }
      }
    }
  }
}`;

/**
 * Injects endpoints for case counts for GENOXA
 */
const genoxaCaseCountSlice = graphqlAPISlice.injectEndpoints({
  endpoints: (builder) => ({
    genoxaCaseCount: builder.query<number, FilterSet>({
      query: (cohortFilters) => {
        const graphQLFilters = {
          cohortFilters: buildCohortGqlOperator(cohortFilters),
          genoxaCaseFilter: GENOXA_FILE_FILTER,
        };
        return {
          graphQLFilters,
          graphQLQuery,
        };
      },
      transformResponse: (response) =>
        response?.data?.viewer?.explore?.genoxaCaseCount?.hits?.total ?? 0,
    }),
  }),
});

export const { useLazyGenoxaCaseCountQuery } = genoxaCaseCountSlice;