---
layout: default
title: Home
nav_order: 1
description: "The core MATLAB Package Manager repository for OceanKit packages"
permalink: /
---

# OceanKit

OceanKit is the core MATLAB Package Manager repository for the OceanKit ecosystem. It is the distribution repository that carries released package snapshots such as `WaveVortexModel-4.0.2`, not the primary authoring location for day-to-day package development.

If you want to install packages, add `OceanKit` as an MPM repository and install the package you need. If you want to author or release a package, work in the package's own Git repository and use the guidance collected here to keep repository structure, MATLAB code, and documentation consistent across the OceanKit ecosystem.

## Start Here

- [Installation](installation) explains how to register the repository and install packages with MPM.
- [Developers guide](developers-guide) collects the package design pattern, MATLAB style guide, and documentation style guide used across OceanKit packages.

## What Lives In OceanKit

- versioned package snapshots for distribution
- reusable release tooling for package repositories
- shared conventions for package authors
