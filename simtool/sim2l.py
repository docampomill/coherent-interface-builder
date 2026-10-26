"""What the Sim2L notebook needs from simtool and yamlmagic, in kernels that have neither (nanoHUB's TorchSim).

With simtool installed, DB and getValidatedInputs are simtool's own. Without it, outputs are written the way
simtool writes them (a scrapbook "text" scrap holding the JSON of the value), so the caller's simtool reads them.
"""
import json

from IPython import get_ipython
from IPython.display import display

try:
    from simtool import DB, getValidatedInputs
except ImportError:
    def getValidatedInputs(inputs):
        """Default values; the caller's simtool has already validated what papermill injects."""
        return {name: spec["value"] for name, spec in inputs.items() if "value" in spec}

    def glue(name, value):
        payload = dict(name=name, data=json.dumps(value), encoder="text", version=1)
        display({"application/scrapbook.scrap.text+json": payload},
                metadata={"scrapbook": dict(name=name, data=True, display=False)}, raw=True)

    class DB:
        def __init__(self, outputs):
            self.outputs, self.unsaved = set(outputs), set(outputs)
            glue("simToolSaveErrorOccurred", 0)
            glue("simToolAllOutputsSaved", 0)

        def save(self, name, value=None, display=False, file=None):
            if name not in self.outputs:
                raise ValueError(f'"{name}" not in output schema!')
            glue(name, value)
            self.unsaved.discard(name)
            if not self.unsaved:
                glue("simToolAllOutputsSaved", 1)


def load_yaml_magic():
    """%%yaml NAME: the cell parsed as YAML into NAME, through yamlmagic if it is installed."""
    shell = get_ipython()
    try:
        shell.run_line_magic("load_ext", "yamlmagic")
    except ModuleNotFoundError:
        import yaml

        def yaml_magic(line, cell):
            shell.user_ns[line.split()[0]] = yaml.safe_load(cell)

        shell.register_magic_function(yaml_magic, magic_kind="cell", magic_name="yaml")
