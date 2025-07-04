import { graphqlAPISlice, GraphQLApiResponse } from "../gdcapi/gdcgraphql";

const geneSetCountQuery = `query geneSetCounts(
  $filters: FiltersArgument
) {
  viewer {
    explore {
      genes {
        hits(filters: $filters, first: 0) {
          total
        }
      }
    }
  }
}
`;

const transformGeneSetCountResponse = (
  response: GraphQLApiResponse<any>,
): number => {
  return response.data.viewer.explore.genes.hits.total;
};

const ssmsSetCountQuery = `query ssmSetCounts(
  $filters: FiltersArgument
) {
  viewer {
    explore {
      ssms {
        hits(filters: $filters, first: 0) {
          total
        }
      }
    }
  }
}`;

const transformSsmsSetCountResponse = (
  response: GraphQLApiResponse<any>,
): number => {
  return response.data.viewer.explore.ssms.hits.total;
};

const caseSetCountQuery = `query caseSetCounts(
  $filters: FiltersArgument
) {
  viewer {
    explore {
      cases {
        hits(filters: $filters, first: 0) {
          total
        }
      }
    }
  }
}`;

const transformCaseSetCountResponse = (
  response: GraphQLApiResponse<any>,
): number => {
  return response.data.viewer.explore.cases.hits.total;
};

export const setCountSlice = graphqlAPISlice
  .enhanceEndpoints({ addTagTypes: ["geneSets", "ssmsSets", "caseSets"] })
  .injectEndpoints({
    endpoints: (builder) => ({
      geneSetCount: builder.query({
        query: ({ setId, additionalFilters }) => ({
          graphQLQuery: geneSetCountQuery,
          graphQLFilters: additionalFilters
            ? {
                filters: {
                  content: [
                    {
                      op: "=",
                      content: {
                        field: "genes.gene_id",
                        value: `set_id:${setId}`,
                      },
                    },
                    additionalFilters,
                  ],
                  op: "and",
                },
              }
            : {
                filters: {
                  op: "=",
                  content: {
                    field: "genes.gene_id",
                    value: `set_id:${setId}`,
                  },
                },
              },
        }),
        transformResponse: transformGeneSetCountResponse,
        providesTags: (_result, _error, arg) => [
          { type: "geneSets", id: arg.setId },
        ],
      }),
      geneSetCounts: builder.query<
        Record<string, number>,
        { setIds: string[] }
      >({
        queryFn: async (_arg, _queryApi, _extraOptions, fetchWithBQ) => {
          const counts: Record<string, number> = {};
          for (const setId of _arg.setIds) {
            const result = await fetchWithBQ({
              graphQLQuery: geneSetCountQuery,
              graphQLFilters: {
                filters: {
                  op: "=",
                  content: {
                    field: "genes.gene_id",
                    value: `set_id:${setId}`,
                  },
                },
              },
            });
            if (result.error) {
              return { error: result.error };
            } else {
              counts[setId as string] = transformGeneSetCountResponse(
                result.data as GraphQLApiResponse<any>,
              );
            }
          }

          return { data: counts };
        },
        providesTags: (_result, _error, arg) =>
          arg.setIds.map((id) => ({ type: "geneSets", id })),
      }),
      ssmSetCount: builder.query({
        query: ({ setId, additionalFilters }) => ({
          graphQLQuery: ssmsSetCountQuery,
          graphQLFilters: additionalFilters
            ? {
                filters: {
                  content: [
                    {
                      op: "=",
                      content: {
                        field: "ssms.ssm_id",
                        value: `set_id:${setId}`,
                      },
                    },
                    additionalFilters,
                  ],
                  op: "and",
                },
              }
            : {
                filters: {
                  op: "=",
                  content: {
                    field: "ssms.ssm_id",
                    value: `set_id:${setId}`,
                  },
                },
              },
        }),
        transformResponse: transformSsmsSetCountResponse,
        providesTags: (_result, _error, arg) => [
          { type: "ssmsSets", id: arg.setId },
        ],
      }),
      ssmSetCounts: builder.query<Record<string, number>, { setIds: string[] }>(
        {
          queryFn: async (_arg, _queryApi, _extraOptions, fetchWithBQ) => {
            const counts: Record<string, number> = {};
            for (const setId of _arg.setIds) {
              const result = await fetchWithBQ({
                graphQLQuery: ssmsSetCountQuery,
                graphQLFilters: {
                  filters: {
                    op: "=",
                    content: {
                      field: "ssms.ssm_id",
                      value: `set_id:${setId}`,
                    },
                  },
                },
              });
              if (result.error) {
                return { error: result.error };
              } else {
                counts[setId as string] = transformSsmsSetCountResponse(
                  result.data as GraphQLApiResponse<any>,
                );
              }
            }

            return { data: counts };
          },
          providesTags: (_result, _error, arg) =>
            arg.setIds.map((id) => ({ type: "ssmsSets", id })),
        },
      ),
      caseSetCount: builder.query({
        query: ({ setId, additionalFilters }) => ({
          graphQLQuery: `query caseSetCounts(
          $filters: FiltersArgument
        ) {
          viewer {
            explore {
              cases {
                hits(filters: $filters, first: 0) {
                  total
                }
              }
            }
          }
        }
        `,
          graphQLFilters: additionalFilters
            ? {
                filters: {
                  content: [
                    {
                      op: "=",
                      content: {
                        field: "cases.case_id",
                        value: `set_id:${setId}`,
                      },
                    },
                    additionalFilters,
                  ],
                  op: "and",
                },
              }
            : {
                filters: {
                  op: "=",
                  content: {
                    field: "cases.case_id",
                    value: `set_id:${setId}`,
                  },
                },
              },
        }),
        transformResponse: transformCaseSetCountResponse,
        providesTags: (_result, _error, arg) => [
          { type: "caseSets", id: arg.setId },
        ],
      }),
      caseSetCounts: builder.query<
        Record<string, number>,
        { setIds: string[] }
      >({
        queryFn: async (_arg, _queryApi, _extraOptions, fetchWithBQ) => {
          const counts: Record<string, number> = {};
          for (const setId of _arg.setIds) {
            const result = await fetchWithBQ({
              graphQLQuery: caseSetCountQuery,
              graphQLFilters: {
                filters: {
                  op: "=",
                  content: {
                    field: "cases.case_id",
                    value: `set_id:${setId}`,
                  },
                },
              },
            });
            if (result.error) {
              return { error: result.error };
            } else {
              counts[setId as string] = transformCaseSetCountResponse(
                result.data as GraphQLApiResponse<any>,
              );
            }
          }
          return { data: counts };
        },
        providesTags: (_result, _error, arg) =>
          arg.setIds.map((id) => ({ type: "caseSets", id })),
      }),
    }),
  });

export const {
  useGeneSetCountQuery,
  useLazyGeneSetCountQuery,
  useGeneSetCountsQuery,
  useSsmSetCountQuery,
  useLazySsmSetCountQuery,
  useSsmSetCountsQuery,
  useCaseSetCountQuery,
  useLazyCaseSetCountQuery,
  useCaseSetCountsQuery,
} = setCountSlice;
