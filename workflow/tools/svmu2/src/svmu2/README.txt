SVMU2 — Trevor D. Millar

SVMU2 is organized using a layered architecture to separate user interface, workflow logic, core algorithms, domain models, and input/output handling. This separation improves reproducibility, testability, and extensibility as the project evolves.

The codebase is divided into five conceptual layers:

Interface
	Entry points that parse user input and translate it into structured program arguments.

Orchestration
	Defines workflow logic. Coordinates models and algorithms to fulfill requests from the interface without implementing core computational logic.

Models
	Domain objects representing key biological and computational concepts (e.g., traversal paths, structural variants, alignment blocks).

Algorithms
	Core computational logic operating on models. Designed to be deterministic, testable, and independent of CLI and file system state.

IO
	Handles input and output operations, including parsing external file formats and writing results. Keeps file-system concerns separate from computational logic.
