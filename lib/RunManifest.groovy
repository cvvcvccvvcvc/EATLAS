import groovy.json.JsonOutput
import groovy.json.JsonSlurper

import java.io.FileOutputStream
import java.nio.charset.StandardCharsets
import java.nio.file.AtomicMoveNotSupportedException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.nio.file.attribute.PosixFilePermissions
import java.security.MessageDigest
import java.util.regex.Pattern


class RunManifest {
    static final int SCHEMA_VERSION = 3
    static final int EVIDENCE_INVENTORY_SCHEMA_VERSION = 1
    static final String REDACTED = '<redacted>'
    private static final Set<String> ACTIVE_RUNS = Collections.synchronizedSet(
        new HashSet<String>()
    )

    private static final Pattern SENSITIVE_KEY = Pattern.compile(
        '(?i).*(password|passwd|secret|token|credential|api[_-]?key|' +
        'access[_-]?key|private[_-]?key).*'
    )
    private static final Pattern URL_USERINFO = Pattern.compile(
        '(?i)([a-z][a-z0-9+.-]*://)[^/@\\s]+@'
    )
    private static final Pattern SENSITIVE_QUERY = Pattern.compile(
        '(?i)([?&](?:password|passwd|secret|token|credential|api[_-]?key|' +
        'access[_-]?key|private[_-]?key)=)[^&#\\s]*'
    )

    static void start(
        Path manifestPath,
        Path projectDir,
        def workflow,
        def parameters,
        Path parameterSchema = null
    ) {
        validateStart(manifestPath, workflow.resume as boolean)
        Map git = gitProvenance(projectDir)
        Map payload = [
            schema_version: SCHEMA_VERSION,
            pipeline: 'gaph_v2',
            status: 'running',
            started_at: workflow.start?.toString(),
            completed_at: null,
            success: null,
            exit_status: null,
            run_name: workflow.runName?.toString(),
            session_id: workflow.sessionId?.toString(),
            profiles: parseProfiles(workflow.profile),
            resume: workflow.resume as boolean,
            nextflow_version: workflow.nextflow?.version?.toString(),
            git_commit: git.commit,
            git_dirty: git.dirty,
            evidence_inventory: null,
            parameters: sanitizeParameters(parameters, parameterSchema),
        ]
        writeAtomic(manifestPath, payload)
        ACTIVE_RUNS.add(runKey(manifestPath, workflow))
    }

    private static void validateStart(Path manifestPath, boolean resume) {
        if (!Files.exists(manifestPath)) {
            return
        }
        if (!Files.isRegularFile(manifestPath)) {
            throw new IllegalStateException(
                "Run manifest is not a regular file: ${manifestPath}"
            )
        }
        Map existing = (Map) new JsonSlurper().parse(manifestPath.toFile())
        if (existing.schema_version != SCHEMA_VERSION) {
            throw new IllegalStateException(
                "Unsupported existing run manifest schema: " +
                "${existing.schema_version} (${manifestPath})"
            )
        }
        if (existing.pipeline != 'gaph_v2') {
            throw new IllegalStateException(
                "Existing output is not a gaph_v2 run: ${manifestPath}"
            )
        }

        String status = existing.status?.toString()
        if (status == 'complete') {
            throw new IllegalStateException(
                "Completed source run is immutable; choose a new --outdir: " +
                manifestPath.parent
            )
        }
        if (!(status in ['running', 'failed'])) {
            throw new IllegalStateException(
                "Existing run manifest has invalid status=${status}: ${manifestPath}"
            )
        }
        if (!resume) {
            throw new IllegalStateException(
                "Existing incomplete run requires -resume: ${manifestPath.parent}"
            )
        }
    }

    static void finish(Path manifestPath, def workflow) {
        if (!ACTIVE_RUNS.remove(runKey(manifestPath, workflow))) {
            return
        }
        if (!Files.isRegularFile(manifestPath)) {
            if (!(workflow.success as boolean)) {
                return
            }
            throw new IllegalStateException("Run manifest not found: ${manifestPath}")
        }
        Map payload = (Map) new JsonSlurper().parse(manifestPath.toFile())
        if (payload.schema_version != SCHEMA_VERSION) {
            throw new IllegalStateException(
                "Unsupported run manifest schema: ${payload.schema_version}"
            )
        }
        String sessionId = workflow.sessionId?.toString()
        if (payload.session_id != sessionId) {
            throw new IllegalStateException(
                "Run manifest session changed: ${payload.session_id} != ${sessionId}"
            )
        }

        boolean success = workflow.success as boolean
        payload.status = success ? 'complete' : 'failed'
        payload.completed_at = workflow.complete?.toString()
        payload.success = success
        payload.exit_status = workflow.exitStatus as Integer
        payload.evidence_inventory = success
            ? evidenceInventoryDescriptor(manifestPath.parent.resolve('evidence_inventory.json'))
            : null
        writeAtomic(manifestPath, payload)
    }

    private static Map evidenceInventoryDescriptor(Path path) {
        if (!Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS)) {
            throw new IllegalStateException("Evidence inventory not found: ${path}")
        }
        byte[] content = Files.readAllBytes(path)
        Map inventory = (Map) new JsonSlurper().parseText(
            new String(content, StandardCharsets.UTF_8)
        )
        if (
            inventory.schema_version != EVIDENCE_INVENTORY_SCHEMA_VERSION ||
            inventory.scope != ['fetch', 'alignment', 'annotation'] ||
            !(inventory.file_count instanceof Number) ||
            (inventory.file_count as Long) < 1L ||
            !(inventory.total_bytes instanceof Number) ||
            (inventory.total_bytes as Long) < 0L ||
            !(inventory.tree_sha256 ==~ /[0-9a-f]{64}/) ||
            !(inventory.files instanceof List)
        ) {
            throw new IllegalStateException("Invalid evidence inventory: ${path}")
        }
        MessageDigest digest = MessageDigest.getInstance('SHA-256')
        return [
            path: 'evidence_inventory.json',
            schema_version: EVIDENCE_INVENTORY_SCHEMA_VERSION,
            size_bytes: content.length,
            sha256: digest.digest(content).encodeHex().toString(),
        ]
    }

    private static String runKey(Path manifestPath, def workflow) {
        return manifestPath.toAbsolutePath().normalize().toString() + '\u0000' +
            workflow.sessionId?.toString()
    }

    private static List<String> parseProfiles(def rawProfiles) {
        String text = rawProfiles?.toString()?.trim()
        if (!text) {
            return []
        }
        return text.split(',')
            .collect { it.trim() }
            .findAll { it }
    }

    private static Map sanitizeParameters(def parameters, Path parameterSchema) {
        Map sanitized = (Map) sanitizeValue(parameters)
        if (parameterSchema == null) {
            return sanitized
        }
        if (!Files.isRegularFile(parameterSchema)) {
            throw new IllegalArgumentException(
                "Parameter schema not found: ${parameterSchema}"
            )
        }
        Map schema = (Map) new JsonSlurper().parse(parameterSchema.toFile())
        Set<String> declared = new HashSet<>()
        Map definitions = (Map) (schema['$defs'] ?: [:])
        definitions.values().each { definition ->
            Map properties = (Map) (definition.properties ?: [:])
            declared.addAll(properties.keySet()*.toString())
        }
        Map<String, Object> selected = new TreeMap<>()
        sanitized.each { key, value ->
            if (declared.contains(key)) {
                selected[key] = value
            }
        }
        return selected
    }

    private static Object sanitizeValue(Object value, String key = null) {
        if (key != null && SENSITIVE_KEY.matcher(key).matches()) {
            return REDACTED
        }
        if (value == null || value instanceof Number || value instanceof Boolean) {
            return value
        }
        if (value instanceof Map) {
            Map<String, Object> sanitized = new TreeMap<>()
            value.each { nestedKey, nestedValue ->
                String name = nestedKey.toString()
                sanitized[name] = sanitizeValue(nestedValue, name)
            }
            return sanitized
        }
        if (value instanceof Collection) {
            return value.collect { sanitizeValue(it) }
        }
        if (value.getClass().isArray()) {
            return Arrays.asList((Object[]) value).collect { sanitizeValue(it) }
        }
        return sanitizeString(value.toString())
    }

    private static String sanitizeString(String value) {
        String withoutUserinfo = URL_USERINFO.matcher(value).replaceAll(
            '$1<redacted>@'
        )
        return SENSITIVE_QUERY.matcher(withoutUserinfo).replaceAll(
            '$1<redacted>'
        )
    }

    private static Map gitProvenance(Path projectDir) {
        try {
            Map commit = runCommand(
                ['git', 'rev-parse', '--verify', 'HEAD'],
                projectDir
            )
            if (commit.exitStatus != 0 || !commit.output) {
                return [commit: null, dirty: null]
            }
            Map status = runCommand(
                ['git', 'status', '--porcelain=v1', '--untracked-files=normal'],
                projectDir
            )
            return [
                commit: commit.output,
                dirty: status.exitStatus == 0 ? !status.output.isEmpty() : null,
            ]
        } catch (IOException ignored) {
            return [commit: null, dirty: null]
        }
    }

    private static Map runCommand(List<String> command, Path directory) {
        Process process = new ProcessBuilder(command)
            .directory(directory.toFile())
            .redirectErrorStream(true)
            .start()
        String output = process.inputStream.getText('UTF-8').trim()
        int exitStatus = process.waitFor()
        return [exitStatus: exitStatus, output: output]
    }

    private static void writeAtomic(Path path, Map payload) {
        Files.createDirectories(path.parent)
        Path temporary = Files.createTempFile(
            path.parent,
            ".${path.fileName}.",
            '.tmp'
        )
        try {
            byte[] content = (
                JsonOutput.prettyPrint(JsonOutput.toJson(payload)) + '\n'
            ).getBytes(StandardCharsets.UTF_8)
            FileOutputStream stream = new FileOutputStream(temporary.toFile())
            try {
                stream.write(content)
                stream.flush()
                stream.fd.sync()
            } finally {
                stream.close()
            }
            try {
                Files.setPosixFilePermissions(
                    temporary,
                    PosixFilePermissions.fromString('rw-r--r--')
                )
            } catch (UnsupportedOperationException ignored) {
                // POSIX permissions are unavailable on some local filesystems.
            }
            try {
                Files.move(
                    temporary,
                    path,
                    StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING
                )
            } catch (AtomicMoveNotSupportedException ignored) {
                Files.move(
                    temporary,
                    path,
                    StandardCopyOption.REPLACE_EXISTING
                )
            }
        } finally {
            Files.deleteIfExists(temporary)
        }
    }
}
