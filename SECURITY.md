# Security Policy

## Supported versions

Covariant / Leonard is currently in alpha. Security fixes are provided on a best-effort basis for the latest published alpha release and the `main` branch.

| Version | Supported |
| --- | --- |
| Latest published alpha | Yes |
| `main` | Yes, best effort |
| Older alpha releases | No |

Because this project is experimental software, users should keep backups of original FlowJo workspace files and should not use Leonard as the sole control for decisions involving sensitive, clinical, regulatory, or production data.

## Reporting a vulnerability

Please do **not** report suspected security vulnerabilities in a public GitHub issue, discussion, pull request, or review comment.

The preferred reporting method is GitHub's private vulnerability reporting feature in the repository's **Security** tab, if it is available.

If private vulnerability reporting is unavailable, contact the maintainer, [@black98fxstc](https://github.com/black98fxstc), privately through GitHub. Please do not include exploit details in a public message. If you need to create an initial public issue to establish contact, provide only a brief statement that a security report requires a private communication channel.

Please include as much of the following information as you can:

- A clear description of the vulnerability and its potential impact.
- The affected release, commit, operating system, and architecture.
- The relevant input type or file, such as a `.wsp` workspace or other data file.
- Reproduction steps or a minimal proof of concept, shared privately.
- Any logs, stack traces, screenshots, or command-line arguments needed to reproduce it.
- Whether the issue involves data disclosure, data modification, arbitrary code execution, denial of service, or another impact.

Please remove personal, patient, proprietary, or otherwise sensitive information from reports and test files whenever possible.

## What to expect

This is a one-person project and response times may vary. As a guideline, the maintainer will try to:

1. Acknowledge a report within 14 days.
2. Investigate and clarify the affected versions and impact.
3. Coordinate a fix or mitigation when practical.
4. Publish an advisory or release note when disclosure is appropriate.

Please allow reasonable time for a fix before making vulnerability details public. Coordinated disclosure helps protect alpha testers and downstream users.

## Scope and known considerations

Leonard is a native desktop application that reads FlowJo workspaces and may rewrite workspace files to add analysis results or gating information. It also launches platform-native file-selection helpers on Windows and macOS. Treat workspaces and downloaded release artifacts as untrusted until they have been reviewed and verified.

Reports are especially useful for issues involving:

- Unsafe handling of workspace-derived filenames, population names, sample names, or other XML content.
- Command or argument injection through input files or file paths.
- Unexpected modification, disclosure, or deletion of user data.
- Unsafe XML, XPath, XSLT, or archive processing.
- Vulnerabilities in the packaged Windows or macOS launchers and installers.
- Dependency vulnerabilities that affect the distributed binaries.

General bugs, crashes, incorrect analysis results, usability issues, and feature requests should be reported through the normal issue tracker instead, unless they also create a security impact.

## Dependency and third-party reports

If a vulnerability is specific to a third-party dependency, please report it to the dependency's maintainers as well as to this project when the vulnerable dependency is bundled into Leonard's release artifacts.

## License and warranty

This security policy does not change the terms of the BSD 3-Clause License. Leonard is alpha software and is provided without warranty; users remain responsible for validating results and protecting their data.
