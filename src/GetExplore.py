from dataclasses import dataclass
import os


@dataclass
class GetExplore:
    expr_file: str

    def __post_init__(self) -> None:
        path = os.path.dirname(self.expr_file)
        base = self.expr_file[:-8]
        export = os.path.join(path, f"{base}-explore.txt")
        with open(export, "w") as file_out:
            file_out.write("[]\n")
            file_out.write("name=\n")

            names = ["expr", "index", "survival", "indexHeader", "info"]
            types = ["expr", "idx", "survival", "ih", "info"]

            for name, type in zip(names, types):
                my_file = f"{base}-{type}.txt"
                filepath = os.path.join(path, my_file)
                file_out.write(f"{name}={filepath}\n")

            file_out.write("key=\n")
            file_out.write(f"source={self.accessionID}")
