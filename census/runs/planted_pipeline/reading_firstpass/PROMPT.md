# Prompt given to the first-pass reader (2026-11-06; post hoc diagnostic)

One Claude agent (claude-opus-5-5) received reading/BRIEF.md and the five input files with this instruction:

"You are one reader in a blinded reading for a methods study. Read the brief, then the input files, and decide each
reply. BRIEF (follow it exactly): census/runs/planted_pipeline/reading/BRIEF.md (It says every reply has the final
verdict KNOWN-VARIANT; here every reply is a first-pass reply whose verdict was KNOWN-VARIANT. Apply the same rule.)
INPUT FILES: reading_firstpass/inputs/C01.md to C05.md. Read only the brief and these five files. Do not open any other
file in the repository, do not search the web, and do not judge whether a reply is right; judge only what each reply
states about the observed behavior on its card. OUTPUT: the JSON list the brief specifies ... to
reading_firstpass/output_blinded.json."

Result: 0 of 13 correct, 2 close (output_blinded.json; KEY.json maps reply ids to run and card).
