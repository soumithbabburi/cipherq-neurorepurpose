# AtomisticSkills Setup Guide

Guide the user step-by-step through setting up AtomisticSkills.

## Step 1: Clone the Repository
Ask the user to run this in their terminal *(Optional: Ask them to fork the repository first if they plan to contribute, and clone their fork instead)*:

```bash
git clone git@github.com:learningmatter-mit/AtomisticSkills.git
cd AtomisticSkills
```

## Step 2: Choose Environments
AtomisticSkills uses separate MCP servers running in different conda environments to manage conflicting MLIP and DFT dependencies.

Present using a list:

**MCP Server Environments (Interactive Tools):**
- [ ] **Base** (`base-agent`): Materials Project queries, VASP I/O, base tools (Highly Recommended)
- [ ] **MACE** (`mace-agent`): MACE models (MP, OMAT, MatPES)
- [ ] **MatGL** (`matgl-agent`): MatGL models (CHGNet, M3GNet) and bandgap prediction
- [ ] **FairChem** (`fairchem-agent`): FairChem models (UMA, ESEN)
- [ ] **Atomate2** (`atomate2-agent`): Remote DFT job management via Jobflow/Jobflow-remote
- [ ] **Smol** (`smol-agent`): Cluster expansion and Monte Carlo
- [ ] **DrugDisc** (`drugdisc-agent`): Drug discovery tools (Fingerprints, Docking, ADMET)
- [ ] **MatterGen** (`mattergen-agent`): Generative crystal design from MatterGen
- [ ] **ADiT** (`adit-agent`): All-atom diffusion generation
- [ ] **DiffCSP** (`diffcsp-agent`): Symmetry-constrained crystal generation

**Script-Only Environments (No MCP Server):**
- [ ] **XRD** (`xrd-agent`): XRD spectrum phase analysis and refinement tools
- [ ] **ORCA** (`orca-agent`): DFT structural optimization and single-points via SCINE wrapper
- [ ] **Phase Field** (`phasefield-agent`): Simulation of grain growth and spinodal decomposition
- [ ] **CALPHAD** (`calphad-agent`): Thermodynamics and phase diagram modeling
- [ ] **NMR** (`nmr-agent`): NMR mixture deconvolution and kinetics prediction
- [ ] **React-OT** (`react-ot-agent`): Transition state structural generation
- [ ] **SCD** (`scd-agent`): Pretrained Self-Conditioned Denoising for property prediction

Options: Keep it simple! Ask exactly which frameworks they want to use. "I only need MACE and basic tools" $\rightarrow$ `base-agent` and `mace-agent`.

## Step 3: Install Conda Environments
Depending on their choices in Step 2, provide the commands to set them up:

```bash
bash conda-envs/base-agent/install.sh
bash conda-envs/mace-agent/install.sh
# ... other selected environments
```
*(Remind the user this might take a few minutes. Wait for them to finish before proceeding).*

## Step 4: Configuration & API Keys
Many tools require API keys (like the Materials Project API) or binary paths.
Have the user create a global configuration file in their home directory.

**Provide this template (`~/.atomistic_skills.yaml`):**
```yaml
# Materials Project API Key (Required for base-server)
MP_API_KEY: "your_mp_api_key_here"

# Atomate2 Remote Project (Required for remote job monitoring)
ATOMATE2_REMOTE_PROJECT: "remote_perlmutter"

# Required for running molecular DFT calculations with ORCA
ORCA_BINARY_PATH: /path/to/orca_directory/orca
```
*(Tell them they can also set these as environment variables like `export MP_API_KEY="key"`, but the file is persistent).*

## Step 5: Configure MCP Servers
The project provides an `mcp_config.json` that defines all tools. The placeholder paths need to be updated to match the user's specific conda path.

Run together:
```bash
python configure_mcp.py
```
*(This auto-detects `miniforge3` or `miniconda3`. If it fails, they can pass the base path manually: `python configure_mcp.py /path/to/miniforge3`)*

## Step 6: Add to AI Assistant
Finally, copy the patched MCP settings into their AI copilot's configuration file.

| Client | Project scope | Global scope (all projects) |
|--------|--------------|----------------------------|
| **Claude Code** | `.mcp.json` (project root) | `~/.claude.json` |
| **Codex CLI** | `.codex/config.toml` | `~/.codex/config.toml` |
| **Cursor** | `.cursor/mcp.json` | `~/.cursor/mcp.json` |
| **Windsurf** | — | `~/.codeium/windsurf/mcp_config.json` |
| **Gemini CLI** | `.gemini/settings.json` | `~/.gemini/settings.json` |

The `configure_mcp.py` script handles path substitution automatically:

```bash
# Project scope only (default) — tools available when inside this repo
python configure_mcp.py

# Global scope — tools available in every project/session
python configure_mcp.py --scope global

# Both scopes
python configure_mcp.py --scope both
```

**Claude Code users: optional global skill awareness**

With `--scope global`, MCP tools become available everywhere. To also make Claude aware of the 129+ skills and workflows when working outside this repo, create `~/.claude/CLAUDE.md`:

```markdown
# AtomisticSkills

When a task involves materials simulation, drug discovery, or atomistic modeling,
the full toolkit lives at /path/to/AtomisticSkills.
Skills: /path/to/AtomisticSkills/.agents/skills/
Workflows: /path/to/AtomisticSkills/.agents/workflows/
Discover skills: grep -r "^description:" /path/to/AtomisticSkills/.agents/skills/*/SKILL.md
```

**Codex CLI users: global MCP, rules, skills, and workflows**

With `--scope global` or `--scope both`, `configure_mcp.py --agent codex` will:
1. Register MCP tools globally in `~/.codex/config.toml`.
2. Append/update the `AtomisticSkills Global Reference` block in `~/.codex/AGENTS.md` for rules and workflows.
3. Create a compact pointer skill at `~/.codex/skills/atomisticskills/SKILL.md`.
4. Symlink every project skill from `.agents/skills/` into `~/.codex/skills/`, making the full skill library available when Codex starts outside this repository.

The symlink step preserves unrelated existing global skills. If a non-project skill with the same name already exists in `~/.codex/skills/`, `configure_mcp.py` skips that entry and reports it. On each run, stale Codex symlinks pointing to removed AtomisticSkills project skills are cleaned up, and new project skills are linked globally. When Codex is started inside the AtomisticSkills repository, the global reference tells it to prefer the project-local `AGENTS.md` and project-local skills, avoiding duplicate rule/skill context.

**Gemini (CLI & IDE/Antigravity) users: global MCP, rules, skills, and workflows**

With `--scope global` or `--scope both`, running `configure_mcp.py --agent gemini` will:
1. Register MCP tools globally in `~/.gemini/config/mcp_config.json` and `~/.gemini/settings.json`.
2. Register a global IDE plugin in `~/.gemini/config/plugins/Google.atomisticskills.atomisticskills` and symlink the skills directory directly, making all 129+ skills natively accessible in any session.
3. Append/update the `AtomisticSkills Global Reference` block in your global rules file (`~/.gemini/GEMINI.md`) for workflows and rules.

This makes Gemini automatically aware of all rules, skills, and workflows when working in any folder outside this repository, while preventing duplicate context/rule loading when working inside the repository. Gemini uses a directory symlink to `.agents/skills/`, so added and removed project skills are reflected globally through that link; rerunning `configure_mcp.py` refreshes the plugin symlink if it becomes stale or points elsewhere.

Restart the assistant after any config changes.

## What's Next? (Guided First Use)
**Run a live test WITH the user.**

**Demo Query (Base agent):** "Search the Materials Project for the stable structure of LiFePO4."
*(Use the `search_materials_project_by_formula` tool)*

If they set up a MLIP (e.g., MACE): "Predict the forces and energy for this LiFePO4 structure using the MACE model."

## Best Practices for Users
- **Leverage Local GPUs**: We highly recommend running the framework on a machine with local GPU resources so MLIP tasks can evaluate quickly without external compute costs.
- **Customize**: Add your own specialized SKILLs, MCP tools, and Workflows directly to the project structure to tailor it to your research needs.
- **Contribute Back**: If you develop a robust, generalized tool or SKILL, please submit a PR to the main branch! We actively acknowledge all open-source contributors.

## Common Issues
| Issue | Fix |
|-------|-----|
| MCP Tools not showing up | Verify JSON syntax in the copilot's config file and restart the IDE/copilot. |
| Tool execution failed / Python not found | Ensure `configure_mcp.py` successfully updated the `command` paths to the correct conda envs. |
| Atomate2 remote worker issues | See `conda-envs/atomate2-agent/atomate2_remote_worker_setup.md` |
| MLIP environment conflicts | Each MCP server handles its own environment isolation automatically via the copied `mcp_config.json`. |
