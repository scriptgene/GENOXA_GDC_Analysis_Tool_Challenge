import { buildCohortGqlOperator, useCurrentCohortFilters } from "@gff/core";
import { GENOXA_FILE_FILTER } from "../genoxa/GenoxaCaseCount";

const SagarPatel_GENOXA = () => {
  const filters = useCurrentCohortFilters();
  const cohortFilter = buildCohortGqlOperator(filters); // undefined if cohort is empty

  // Always include GENOXA_FILE_FILTER. Only add the cohort filter if one exists.
  const combinedFilter = cohortFilter
    ? { op: "and", content: [GENOXA_FILE_FILTER, cohortFilter] }
    : GENOXA_FILE_FILTER;

  const src = `http://localhost:8081/?filters=${encodeURIComponent(
    JSON.stringify(combinedFilter),
  )}`;

  return (
    <div style={{ height: "100vh" }}>
      <iframe
        key={src}
        src={src}
        style={{ width: "100%", height: "100%", border: "none" }}
        title="Genoxa"
      />
    </div>
  );
};

export default SagarPatel_GENOXA;