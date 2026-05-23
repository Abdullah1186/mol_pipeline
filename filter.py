"""Back-compat shim. Real pipeline lives in pipelines/. Edit pipelines/params.yaml to configure."""
from pipelines.runner import main

if __name__ == "__main__":
    main()
