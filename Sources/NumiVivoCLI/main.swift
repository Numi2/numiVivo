import Foundation
import NumiVivoKit

private let arguments = Array(CommandLine.arguments.dropFirst())

switch arguments.first {
case "--version", "version":
    print("NumiVivo \(NumiVivo.version) (program pack v\(NumiVivo.programPackVersion))")
case "--help", "help", nil:
    print("""
    NumiVivo \(NumiVivo.version)

    Usage:
      numivivo compile <program.json> --output <program.vivopack>
      numivivo verify <program.json>
      numivivo inspect <program.vivopack>
      numivivo simulate <program.vivopack> --experiment <experiment.json>
      numivivo export-sbml <program.json> --output <model.xml>
      numivivo export-sbol <program.json> --output <design.jsonld>
      numivivo version

    The foundational command surface is installed. Command implementations are
    provided by NumiVivoKit as the compiler and runtime modules are enabled.
    """)
default:
    fputs("numivivo: unknown command '\(arguments[0])'\nRun 'numivivo help' for usage.\n", stderr)
    exit(64)
}
