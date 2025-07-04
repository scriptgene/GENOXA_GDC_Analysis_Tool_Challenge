import React, { ComponentType, ReactElement } from "react";
import { render, RenderOptions, RenderResult } from "@testing-library/react";
import { CoreProvider } from "@gff/core";
import { createTheme, MantineProvider } from "@mantine/core";
import { SummaryModalContext, URLContext } from "src/utils/contexts";
import { DashboardDownloadContext } from "@gff/portal-components";
import tailwindConfig from "tailwind.config";

const AllTheProviders = ({ children }: { children: React.ReactNode }) => {
  const theme = createTheme(
    Object.fromEntries(
      Object.entries(
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        tailwindConfig.plugins.slice(-1)[0].__options.defaultTheme.extend
          .colors,
      ).map(([key, values]) => [key, Object.values(values)]),
    ),
  );

  return (
    <CoreProvider>
      <MantineProvider theme={theme}>
        <URLContext.Provider value={{ prevPath: "", currentPath: "" }}>
          <SummaryModalContext.Provider
            value={{
              entityMetadata: {
                entity_type: null,
                entity_id: null,
              },
              setEntityMetadata: jest.fn(),
            }}
          >
            <DashboardDownloadContext.Provider
              value={{ state: [], dispatch: jest.fn() }}
            >
              {children}
            </DashboardDownloadContext.Provider>
          </SummaryModalContext.Provider>
        </URLContext.Provider>
      </MantineProvider>
    </CoreProvider>
  );
};

const customRender = (
  ui: ReactElement,
  options?: Omit<RenderOptions, "wrapper">,
): RenderResult =>
  render(ui, { wrapper: AllTheProviders as ComponentType, ...options });

export * from "@testing-library/react";
export { customRender as render };
