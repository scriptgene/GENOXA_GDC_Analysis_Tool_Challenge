import userEvent from "@testing-library/user-event";
import * as facetHooks from "../../facets/hooks";
import * as router from "next/router";
import { render } from "test-utils";
import CDaveCard from "./CDaveCard";
import {
  selectFacetDefinitionByName,
  useGetContinuousDataStatsQuery,
} from "@gff/core";

jest.mock("@gff/core", () => ({
  ...jest.requireActual("@gff/core"),
  selectFacetDefinitionByName: jest.fn(),
  useGetContinuousDataStatsQuery: jest.fn(),
}));

jest.spyOn(router, "useRouter").mockImplementation(
  () =>
    ({
      pathname: "",
      query: {},
    } as any),
);

describe("CDaveCard", () => {
  it("enum result with data", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "demographic.gender",
          type: "keyword",
        } as any),
    );
    const data = {
      buckets: [
        { doc_count: 1000, key: "male" },
        { doc_count: 500, key: "female" },
      ],
    };

    const { getByRole } = render(
      <CDaveCard
        data={data}
        field={"demographic.gender"}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(
      getByRole("row", {
        name: "Select male male 1,000 (66.67%)",
      }),
    ).toBeInTheDocument();
  });

  it("enum result with only missing values", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "demographic.gender",
          type: "keyword",
        } as any),
    );
    const data = {
      buckets: [{ doc_count: 1000, key: "_missing" }],
    };

    const { getByText } = render(
      <CDaveCard
        data={data}
        field={"demographic.gender"}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(getByText("No data for this property")).toBeInTheDocument();
  });

  it("enum results with a missing value", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "demographic.gender",
          type: "keyword",
        } as any),
    );
    const data = {
      buckets: [
        { doc_count: 1000, key: "_missing" },
        { doc_count: 500, key: "female" },
      ],
    };

    const { getByRole } = render(
      <CDaveCard
        data={data}
        field={"demographic.gender"}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(
      getByRole("row", {
        name: "Select missing missing 1,000 (66.67%)",
      }),
    ).toBeInTheDocument();
  });

  it("categorical results sorted by count", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "demographic.gender",
          type: "keyword",
        } as any),
    );
    const data = {
      buckets: [
        { doc_count: 500, key: "female" },
        { doc_count: 1000, key: "_missing" },
      ],
    };

    const { getAllByRole } = render(
      <CDaveCard
        data={data}
        field={"demographic.gender"}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(getAllByRole("row")[1]).toHaveTextContent("missing");
    expect(getAllByRole("row")[2]).toHaveTextContent("female");
  });

  it("continuous result with data", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "exposures.cigarettes_per_day",
          type: "long",
        } as any),
    );
    jest.spyOn(facetHooks, "useRangeFacet").mockReturnValue({
      data: { "0.0-12.0": 10, "12.0-24.0": 90 },
      isFetching: false,
      isSuccess: true,
    } as any);
    jest.mocked(useGetContinuousDataStatsQuery).mockReturnValue({
      isFetching: false,
      isSuccess: true,
    } as any);

    const { getByRole } = render(
      <CDaveCard
        data={{ stats: { count: 100 } } as any}
        field={"exposures.cigarettes_per_day"}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(
      getByRole("row", {
        name: "Select 0 to <12 0 to <12 10 (10.00%)",
      }),
    ).toBeInTheDocument();
    expect(
      getByRole("row", {
        name: "Select 12 to <24 12 to <24 90 (90.00%)",
      }),
    ).toBeInTheDocument();
  });

  it("continuous result with negative bucket", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "exposures.cigarettes_per_day",
          type: "long",
        } as any),
    );
    jest.spyOn(facetHooks, "useRangeFacet").mockReturnValue({
      data: { "-28.0-166.8000001": 38 },
      isFetching: false,
      isSuccess: true,
    } as any);
    jest.mocked(useGetContinuousDataStatsQuery).mockReturnValue({
      isFetching: false,
      isSuccess: true,
    } as any);

    const { getByRole } = render(
      <CDaveCard
        data={{ stats: { count: 38 } } as any}
        field="exposures.cigarettes_per_day"
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(getByRole("cell", { name: "-28 to <167" })).toBeInTheDocument();
  });

  it("continuous result with toggled value bucket", async () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "diagnoses.treatments.days_to_treatment_start",
          type: "long",
        } as any),
    );
    jest.spyOn(facetHooks, "useRangeFacet").mockReturnValue({
      data: { "7201.0-12255.8": 10, "12255.8-17310.6": 90 },
      isFetching: false,
      isSuccess: true,
    } as any);
    jest.mocked(useGetContinuousDataStatsQuery).mockReturnValue({
      isFetching: false,
      isSuccess: true,
    } as any);
    jest.spyOn(router, "useRouter").mockImplementation(
      () =>
        ({
          pathname: "",
          query: { featureFlag: "yearToggle" },
        } as any),
    );

    const { getByRole, getByLabelText } = render(
      <CDaveCard
        data={{ stats: { count: 38 } } as any}
        field="diagnoses.treatments.days_to_treatment_start"
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(getByRole("cell", { name: "20 to <34" })).toBeInTheDocument();
    expect(getByRole("cell", { name: "34 to <47" })).toBeInTheDocument();

    await userEvent.click(getByLabelText("Days"));

    expect(getByRole("cell", { name: "7,201 to <12,256" })).toBeInTheDocument();
    expect(
      getByRole("cell", { name: "12,256 to <17,311" }),
    ).toBeInTheDocument();
  });

  it("continuous result with no data", () => {
    jest.mocked(selectFacetDefinitionByName).mockImplementation(
      () =>
        ({
          field: "diagnoses.treatments.days_to_treatment_start",
          type: "long",
        } as any),
    );
    jest
      .spyOn(facetHooks, "useRangeFacet")
      .mockReturnValue({ data: {}, isFetching: false } as any);

    const stats = {
      stats: {
        count: 0,
        min: null,
        max: null,
        avg: null,
        sum: null,
      },
    };

    const { getByText } = render(
      <CDaveCard
        data={stats}
        field={""}
        updateFields={jest.fn()}
        initialDashboardRender
        cohortFilters={undefined}
      />,
    );

    expect(getByText("No data for this property")).toBeInTheDocument();
  });
});
