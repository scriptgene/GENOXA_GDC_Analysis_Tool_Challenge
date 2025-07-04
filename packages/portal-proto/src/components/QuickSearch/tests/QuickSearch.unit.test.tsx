import { waitFor } from "@testing-library/react";
import { render } from "test-utils";
import userEvent from "@testing-library/user-event";
import * as router from "next/router";
import { QuickSearch } from "../QuickSearch";
import { useGetHistoryQuery, useQuickSearchQuery } from "@gff/core";

jest.mock("@gff/core", () => ({
  ...jest.requireActual("@gff/core"),
  useQuickSearchQuery: jest.fn(),
  useGetHistoryQuery: jest.fn(),
}));

jest.spyOn(router, "useRouter").mockImplementation(
  () =>
    ({
      pathname: "",
      query: {},
    } as any),
);

describe("<QuickSearch />", () => {
  test("displays results", async () => {
    jest.mocked(useQuickSearchQuery).mockReturnValue({
      data: {
        searchList: [
          { id: btoa("Gene:111"), symbol: "TH" },
          { id: btoa("Project:222"), name: "Thymoma" },
        ],
        query: "th",
      },
    } as any);
    jest.mocked(useGetHistoryQuery).mockReturnValue({
      data: undefined,
      isSuccess: true,
      isUninitialized: false,
    } as any);

    const { getAllByTestId, getByTestId } = render(<QuickSearch />);
    userEvent.click(getByTestId("textbox-quick-search-bar"));
    userEvent.type(getByTestId("textbox-quick-search-bar"), "th");

    await waitFor(
      () =>
        expect(getAllByTestId("text-search-result").length).toBeGreaterThan(0),
      { timeout: 3000 },
    );
    const results = getAllByTestId("text-search-result");
    expect(results[0].textContent).toEqual("THTHCategory: Gene");
    expect(results[1].textContent).toEqual("222ThymomaCategory: Project");
  });

  test("displays superseded file", async () => {
    jest.mocked(useQuickSearchQuery).mockImplementation(
      (search_str) =>
        (search_str === "111-222"
          ? {
              data: {
                searchList: [],
                query: "111-222",
              },
            }
          : {
              data: {
                searchList: [{ uuid: "444-555" }],
                query: "444-555",
              },
            }) as any,
    );
    jest.mocked(useGetHistoryQuery).mockReturnValue({
      data: [
        {
          uuid: "444-555",
          file_change: "released",
        },
        { uuid: "111-222", file_change: "superseded" },
      ],
      isSuccess: true,
      isUninitialized: false,
    } as any);

    const { getAllByTestId, getByTestId } = render(<QuickSearch />);
    userEvent.click(getByTestId("textbox-quick-search-bar"));
    userEvent.type(getByTestId("textbox-quick-search-bar"), "111-222");

    await waitFor(
      () =>
        expect(getAllByTestId("text-search-result").length).toBeGreaterThan(0),
      { timeout: 3000 },
    );
    const results = getAllByTestId("text-search-result");
    expect(results[0].textContent).toEqual(
      "444-555File 111-222 has been updatedCategory: File",
    );
  });

  test("displays no results found if history is empty", async () => {
    jest.mocked(useQuickSearchQuery).mockReturnValue({
      data: {
        searchList: [],
        query: "111-222",
      },
    } as any);
    jest.mocked(useGetHistoryQuery).mockReturnValue({
      data: [],
      isSuccess: true,
      isUninitialized: false,
    } as any);

    const { getByTestId } = render(<QuickSearch />);
    userEvent.click(getByTestId("textbox-quick-search-bar"));
    userEvent.type(getByTestId("textbox-quick-search-bar"), "111-222");

    await waitFor(
      () =>
        expect(getByTestId("no-results-quick-search-bar")).toBeInTheDocument(),
      { timeout: 3000 },
    );
    const noResults = getByTestId("no-results-quick-search-bar");
    expect(noResults).toBeInTheDocument();
  });

  test("displays no results found if file doesn't have a newer version", async () => {
    jest.mocked(useQuickSearchQuery).mockReturnValue({
      data: {
        searchList: [],
        query: "111-222",
      },
    } as any);
    jest.mocked(useGetHistoryQuery).mockReturnValue({
      data: [
        {
          file_change: "released",
          uuid: "111-222",
          version: "1",
        },
      ],
      isSuccess: true,
      isUninitialized: false,
    } as any);

    const { getByTestId } = render(<QuickSearch />);
    userEvent.click(getByTestId("textbox-quick-search-bar"));
    userEvent.type(getByTestId("textbox-quick-search-bar"), "111-222");

    await waitFor(
      () =>
        expect(getByTestId("no-results-quick-search-bar")).toBeInTheDocument(),
      { timeout: 3000 },
    );
    const noResults = getByTestId("no-results-quick-search-bar");
    expect(noResults).toBeInTheDocument();
  });

  test("displays no results in superseded file does not exist", async () => {
    jest.mocked(useQuickSearchQuery).mockImplementation(
      (search_str) =>
        (search_str === "111-222"
          ? {
              data: {
                searchList: [],
                query: "111-222",
              },
            }
          : {
              data: {
                searchList: [],
                query: "444-555",
              },
              isSuccess: true,
            }) as any,
    );
    jest.mocked(useGetHistoryQuery).mockReturnValue({
      data: [
        {
          uuid: "444-555",
          file_change: "released",
        },
        { uuid: "111-222", file_change: "superseded" },
      ],
      isSuccess: true,
      isUninitialized: false,
    } as any);

    const { getByTestId } = render(<QuickSearch />);
    userEvent.click(getByTestId("textbox-quick-search-bar"));
    userEvent.type(getByTestId("textbox-quick-search-bar"), "111-222");

    await waitFor(
      () =>
        expect(getByTestId("no-results-quick-search-bar")).toBeInTheDocument(),
      { timeout: 5000 },
    );
    const noResults = getByTestId("no-results-quick-search-bar");
    expect(noResults).toBeInTheDocument();
  });
});
